# tRNAscanSE/CM.pm
# This class contains parameters and functions for running CM tRNA search used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::CM;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::ScanResult;
use tRNAscanSE::SS;
use tRNAscanSE::Sequence;

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
    
    $self->{CM_mode} = "cove";
    
    $self->{cm_cutoff} = 20;                    # default cutoff score for covels reporting of tRNA
    
    $self->{nci_scan_cutoff} = 70;              # default cutoff score for rescanning noncanonical introns
    
    $self->{split_tRNA_scan_cutoff} = 38;       # default cutoff score for rescanning split tRNA

    $self->{half_tRNA_cutoff} = 15;             # default cutoff score for half tRNA
    
    $self->{BHB_cm_cutoff} = 5.5;               # default score for considering non-canonical intron
    
    $self->{max_tRNA_length} = 500;             # max size of -w parameter passed to covels
                                                #  when using a pre-scanner (eufind or tRNAscan)
                                                
    $self->{max_cove_tRNA_length} = 250;        # max size of -w param if only 
                                                #  Cove is being used (too slow otherwise)
                                                
    $self->{max_cmsearch_tRNA_length} = 250;    # max size of -w param if only 
                                                #  cmsearch is being used (too slow otherwise)

    $self->{CM_check_for_introns} = 0;          # check for non-canonical introns

    $self->{CM_check_for_split_halves} = 0;     # check for split tRNA fragments
                                                
    $self->{min_tRNA_no_intron} = 76;           # min length for average tRNA with no intron
    
    $self->{left_splicing_len} = 27;
    $self->{right_splicing_len} = 28;
    
    $self->{min_intron_length} = 5;             # min size of introns detected by parsing of 
                                                #  coves output

    $self->{skip_pseudo_filter} = 0;            # enable filter for psuedogenes (Cove score <40,
                                                # primary struct score <10 bits, secondary 
                                                # structure score < 5 bits)
                                                
    $self->{min_cove_pseudo_filter_score} = 55; # Below this score, tRNAs are checked
                                                # for min primary and secondary structure
                                                # scores to catch pseudogene repeats
                                                # like rat ID & rodent B2 elements
                                                
     $self->{min_cmsearch_pseudo_filter_score} = 55;     # Below this score, tRNAs are checked
                                                         # for min primary and secondary structure
                                                         # scores to catch pseudogene repeats
                                                         # like rat ID & rodent B2 elements
                                               
    $self->{min_ss_score} = 5;                  # Below this secondary structure score,
                                                #  tRNA is considered a pseudogene
                                                
    $self->{min_hmm_score} = 10;                # Below this primary structure score,
                                                #  tRNA is considered a pseudogene

    $self->{get_hmm_score} = 0;                 # also score tRNA with covariance model
                                                #  without sec structure info, similar
                                                #  to getting hmm score for match of 
                                                #  seq to tRNA hmm  (-H option)
                                                
    $self->{alt_cm_file} = '';                  # alternate covariance model file (-c option)
    
    $self->{main_cm_file} = '';                 # Convariance model file name
    $self->{mainNS_cm_file} = '';
    $self->{arch_gw_scan_cm_file} = '';
    $self->{arch_intron_cm_file} = '';
    $self->{arch_five_half_cm_file} = '';
    $self->{arch_three_half_cm_file} = '';
    $self->{Pselc_cm_file} = '';
    $self->{Eselc_cm_file} = '';
    
    $self->{main_cm_file_path} = '';            # Convariance model file path
    $self->{mainNS_cm_file_path} = '';
    $self->{arch_gw_scan_cm_file_path} = '';
    $self->{arch_intron_cm_file_path} = '';
    $self->{arch_five_half_cm_file_path} = '';
    $self->{arch_three_half_cm_file_path} = '';
    $self->{Pselc_cm_file_path} = '';
    $self->{Eselc_cm_file_path} = '';

    $self->{covels_bin} = "covels-SE";          # Application executable name
    $self->{coves_bin} = "coves-SE";
    $self->{cmsearch_bin} = "cmsearch";

    $self->{tab_results} = +[];
}

sub CM_mode
{
    my $self = shift;
    if (@_) { $self->{CM_mode} = shift; }
    return $self->{CM_mode};
}

sub cove_mode
{
    my $self = shift;
    return ($self->{CM_mode} eq 'cove');
}

sub infernal_mode
{
    my $self = shift;
    return ($self->{CM_mode} eq 'infernal');
}

sub cm_cutoff
{
    my $self = shift;
    if (@_) { $self->{cm_cutoff} = shift; }
    return $self->{cm_cutoff};
}

sub BHB_cm_cutoff
{
    my $self = shift;
    if (@_) { $self->{BHB_cm_cutoff} = shift; }
    return $self->{BHB_cm_cutoff};
}

sub max_tRNA_length
{
    my $self = shift;
    if (@_) { $self->{max_tRNA_length} = shift; }
    return $self->{max_tRNA_length};
}

sub max_cove_tRNA_length
{
    my $self = shift;
    if (@_) { $self->{max_cove_tRNA_length} = shift; }
    return $self->{max_cove_tRNA_length};
}

sub max_cmsearch_tRNA_length
{
    my $self = shift;
    if (@_) { $self->{max_cmsearch_tRNA_length} = shift; }
    return $self->{max_cmsearch_tRNA_length};
}

sub CM_check_for_introns
{
    my $self = shift;
    if (@_) { $self->{CM_check_for_introns} = shift; }
    return $self->{CM_check_for_introns};
}

sub CM_check_for_split_halves
{
    my $self = shift;
    if (@_) { $self->{CM_check_for_split_halves} = shift; }
    return $self->{CM_check_for_split_halves};
}

sub min_tRNA_no_intron
{
    my $self = shift;
    if (@_) { $self->{min_tRNA_no_intron} = shift; }
    return $self->{min_tRNA_no_intron};
}

sub min_intron_length
{
    my $self = shift;
    if (@_) { $self->{min_intron_length} = shift; }
    return $self->{min_intron_length};
}

sub skip_pseudo_filter
{
    my $self = shift;
    if (@_) { $self->{skip_pseudo_filter} = shift; }
    return $self->{skip_pseudo_filter};
}

sub min_pseudo_filter_score
{
    my $self = shift;
    if (@_) { $self->{min_pseudo_filter_score} = shift; }
    return $self->{min_pseudo_filter_score};
}

sub min_ss_score
{
    my $self = shift;
    if (@_) { $self->{min_ss_score} = shift; }
    return $self->{min_ss_score};
}

sub min_hmm_score
{
    my $self = shift;
    if (@_) { $self->{min_hmm_score} = shift; }
    return $self->{min_hmm_score};
}

sub get_hmm_score
{
    my $self = shift;
    if (@_) { $self->{get_hmm_score} = shift; }
    return $self->{get_hmm_score};
}

sub alt_cm_file
{
    my $self = shift;
    if (@_) { $self->{alt_cm_file} = shift; }
    return $self->{alt_cm_file};
}

sub main_cm_file
{
    my $self = shift;
    if (@_) { $self->{main_cm_file} = shift; }
    return $self->{main_cm_file};
}

sub mainNS_cm_file
{
    my $self = shift;
    if (@_) { $self->{mainNS_cm_file} = shift; }
    return $self->{mainNS_cm_file};
}

sub arch_intron_cm_file
{
    my $self = shift;
    if (@_) { $self->{arch_intron_cm_file} = shift; }
    return $self->{arch_intron_cm_file};
}

sub Pselc_cm_file
{
    my $self = shift;
    if (@_) { $self->{Pselc_cm_file} = shift; }
    return $self->{Pselc_cm_file};
}

sub Eselc_cm_file
{
    my $self = shift;
    if (@_) { $self->{Eselc_cm_file} = shift; }
    return $self->{Eselc_cm_file};
}

sub main_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{main_cm_file_path} = shift; }
    return $self->{main_cm_file_path};
}

sub mainNS_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{mainNS_cm_file_path} = shift; }
    return $self->{mainNS_cm_file_path};
}

sub arch_intron_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{arch_intron_cm_file_path} = shift; }
    return $self->{arch_intron_cm_file_path};
}

sub Pselc_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{Pselc_cm_file_path} = shift; }
    return $self->{Pselc_cm_file_path};
}

sub Eselc_cm_file_path
{
    my $self = shift;
    if (@_) { $self->{Eselc_cm_file_path} = shift; }
    return $self->{Eselc_cm_file_path};
}

sub covels_bin
{
    my $self = shift;
    if (@_) { $self->{covels_bin} = shift; }
    return $self->{covels_bin};
}

sub coves_bin
{
    my $self = shift;
    if (@_) { $self->{coves_bin} = shift; }
    return $self->{coves_bin};
}

sub cmsearch_bin
{
    my $self = shift;
    if (@_) { $self->{cmsearch_bin} = shift; }
    return $self->{cmsearch_bin};
}

sub tab_results
{
    my $self = shift;
    if (@_) { $self->{tab_results} = shift; }
    return $self->{tab_results};
}

sub set_file_paths {
    
    my $self = shift;
    my $opts = shift;
    
    if ($opts->general_mode()) {
        if ($self->infernal_mode()) {
            $self->{main_cm_file} =   "TRNAinf-c.cm";                # use original covariance model 
            $self->{mainNS_cm_file} = "TRNAinf-ns-c.cm";             # no sec struct
        }
        elsif ($self->cove_mode()) {
            $self->{main_cm_file} =   "TRNA2.cm";                   # use original covariance model 
            $self->{mainNS_cm_file} = "TRNA2ns.cm";                 # no sec struct
        }
    }
    elsif ($opts->bact_mode()) {
        if ($self->infernal_mode()) {
            $self->{main_cm_file} =   "TRNAinf-bact-c.cm";           # use bacterial covariance model 
            $self->{mainNS_cm_file} = "TRNAinf-bact-ns-c.cm";        # no sec struct
        }
        elsif ($self->cove_mode()) {
            $self->{main_cm_file} =   "TRNA2-bact.cm";              # use bacterial covariance model 
            $self->{mainNS_cm_file} = "TRNA2-bactns.cm";            # no sec struct
        }
    }
    elsif ($opts->arch_mode()) {
        $self->{arch_intron_cm_file} = "Archaea-BHB-noncan.cm";     # model for finding noncanonical tRNAs
        $self->{arch_five_half_cm_file} = "TRNAinf-arch-5h-nc.cm";         # model for finding 5'half
        $self->{arch_three_half_cm_file} = "TRNAinf-arch-3h-nc.cm";         # model for finding 3'half
        $self->{arch_gw_scan_cm_file} = 'TRNAinf-arch-c.cm';
        if ($self->infernal_mode()) {
            $self->{main_cm_file} =   "TRNAinf-arch-c.cm";           # use archae covariance model 
            $self->{mainNS_cm_file} = "TRNAinf-arch-ns-c.cm";        # no sec struct
        }
        elsif ($opts->cove_mode()) {
            $self->{main_cm_file} =   "TRNA2-arch.cm";              # use archae covariance model 
            $self->{mainNS_cm_file} = "TRNA2-archns.cm";            # no sec struct
        }
    }
    else {
        if ($self->infernal_mode()) {
            $self->{main_cm_file} =   "TRNAinf-euk-c.cm";            # default to eukar cove model 
            $self->{mainNS_cm_file} = "TRNAinf-euk-ns-c.cm";         # no secondary struct
        }
        elsif ($self->cove_mode()) {
            $self->{main_cm_file} =   "TRNA2-euk.cm";               # default to eukar cove model 
            $self->{mainNS_cm_file} = "TRNA2-eukns.cm";             # no secondary struct
        }
    }                           
    
    if ($self->{alt_cm_file} ne '') {
        $self->{main_cm_file} = $self->{Alt_cm_file};               # use alternate cm file specified
                                                                    #  on command line with -c param
        if ($self->infernal_mode()) {
            $self->{mainNS_cm_file} = "TRNAinf-ns-c.cm";
        }
        elsif ($self->cove_mode()) {
            $self->{mainNS_cm_file} = "TRNA2ns.cm";
        }
    }
        
    if ($self->infernal_mode()) {
        $self->{Pselc_cm_file} = "PSELCinf-c.cm";
        $self->{Eselc_cm_file} = "ESELCinf-c.cm";
    }
    elsif ($self->cove_mode()) {
        $self->{Pselc_cm_file} = "PSELC.cm";
        $self->{Eselc_cm_file} = "ESELC.cm";
    }
}

sub check_lib_files {
    
    my $self = shift;
    my $opts = shift;
    my $lib_dir = shift;
    
    if (-r $self->{main_cm_file}) {
        $self->{main_cm_file_path} = $self->{main_cm_file};
    }
    elsif (-r $lib_dir.$self->{main_cm_file}) {
        $self->{main_cm_file_path} =  $lib_dir.$self->{main_cm_file}; 
    }
    else {
        die "FATAL: Unable to open ".$self->{main_cm_file}." covariance model file\n\n";
    }

    if (-r $self->{mainNS_cm_file}) {
        $self->{mainNS_cm_file_path} = $self->{mainNS_cm_file};
    }
    elsif (-r $lib_dir.$self->{mainNS_cm_file}) {
        $self->{mainNS_cm_file_path} =  $lib_dir.$self->{mainNS_cm_file}; 
    }
    else {
        die "FATAL: Unable to open ".$self->{mainNS_cm_file}." covariance model file\n\n";
    }

    if (-r $self->{Pselc_cm_file}) {
        $self->{Pselc_cm_file_path} = $self->{Pselc_cm_file};
    }
    elsif (-r  $lib_dir.$self->{Pselc_cm_file}) {
        $self->{Pselc_cm_file_path} =  $lib_dir.$self->{Pselc_cm_file}; 
    }
    else {
        die "FATAL: Unable to open ".$self->{Pselc_cm_file}." covariance model file\n\n";
    }

    if (-r $self->{Eselc_cm_file}) {
        $self->{Eselc_cm_file_path} = $self->{Eselc_cm_file};
    }
    elsif (-r  $lib_dir.$self->{Eselc_cm_file}) {
        $self->{Eselc_cm_file_path} =  $lib_dir.$self->{Eselc_cm_file}; 
    }
    else {
        die "FATAL: Unable to open ".$self->{Eselc_cm_file}." covariance model file\n\n";
    }
    if ($opts->arch_mode() && ($self->infernal_mode() ||  $self->cove_mode())) {
        if (-r $self->{arch_gw_scan_cm_file}) {
            $self->{arch_gw_scan_cm_file_path} = $self->{arch_gw_scan_cm_file};
        }
        elsif (-r  $lib_dir.$self->{arch_gw_scan_cm_file}) {
            $self->{arch_gw_scan_cm_file_path} =  $lib_dir.$self->{arch_gw_scan_cm_file}; 
        }
        else {
            die "FATAL: Unable to open ".$self->{arch_gw_scan_cm_file}." covariance model file\n\n";
        }
        if (-r $self->{arch_intron_cm_file}) {
            $self->{arch_intron_cm_file_path} = $self->{arch_intron_cm_file};
        }
        elsif (-r  $lib_dir.$self->{arch_intron_cm_file}) {
            $self->{arch_intron_cm_file_path} =  $lib_dir.$self->{arch_intron_cm_file}; 
        }
        else {
            die "FATAL: Unable to open ".$self->{arch_intron_cm_file}." covariance model file\n\n";
        }
        if (-r $self->{arch_five_half_cm_file}) {
            $self->{arch_five_half_cm_file_path} = $self->{arch_five_half_cm_file};
        }
        elsif (-r  $lib_dir.$self->{arch_five_half_cm_file}) {
            $self->{arch_five_half_cm_file_path} =  $lib_dir.$self->{arch_five_half_cm_file}; 
        }
        else {
            die "FATAL: Unable to open ".$self->{arch_five_half_cm_file}." covariance model file\n\n";
        }
        if (-r $self->{arch_three_half_cm_file}) {
            $self->{arch_three_half_cm_file_path} = $self->{arch_three_half_cm_file};
        }
        elsif (-r  $lib_dir.$self->{arch_three_half_cm_file}) {
            $self->{arch_three_half_cm_file_path} =  $lib_dir.$self->{arch_three_half_cm_file}; 
        }
        else {
            die "FATAL: Unable to open ".$self->{arch_three_half_cm_file}." covariance model file\n\n";
        }
     }
}

sub set_bin {
    
    my $self = shift;
    my $bindir = shift;
    
    if ($^O =~ /^MSWin/) {
        $self->{cmsearch_bin} .= ".exe";
        $self->{covels_bin} .= ".exe";
        $self->{coves_bin} .= ".exe";
    }
    if ($self->infernal_mode()) {
        if (!(-x $self->{cmsearch_bin})) {
            $self->{cmsearch_bin} = $bindir.$self->{cmsearch_bin};
            if (!(-x $self->{cmsearch_bin})) {
                die "FATAL: Unable to find ".$self->{cmsearch_bin}." executable\n\n";
            }
        }
    }
    if ($self->cove_mode()) {
        if (!(-x $self->{covels_bin})) {
            $self->{covels_bin} = $bindir.$self->{covels_bin};
            if (!(-x $self->{covels_bin})) {
                die "FATAL: Unable to find ".$self->{covels_bin}." executable\n\n";
            }
        }
        if (!(-x $self->{coves_bin})) {
            $self->{coves_bin} = $bindir.$self->{coves_bin};
            if (!(-x $self->{coves_bin})) {
                die "FATAL: Unable to find ".$self->{coves_bin}." executable\n\n";
            }
        }
    }
}

sub set_search_params {

    my $self = shift;
    my $opts = shift;
    my ($r_scan_len, $r_cur_cm_file,
           $max_search_tRNA_length, $trna_len, $trna_isotype, $ns_cm) = @_;

    # don't set '-W' param over 200 bp if a pre-scanner is being used,
    #  use max window of 150 bp if cmsearch only (too slow otherwise)
    
    if ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run()) {
        $$r_scan_len = &min($trna_len, $self->{max_tRNA_length});
    }
    else {
        $$r_scan_len = $max_search_tRNA_length;
    }        
    
    # set correct CM file for current tRNA
    if ($ns_cm) {
        $$r_cur_cm_file = $self->{mainNS_cm_file_path};
    }
    else {
        $$r_cur_cm_file = $self->{main_cm_file_path};
    
        if ($opts->eufind_mode()) {
            if ($trna_isotype eq "SeCp") {       # use arch/prok selcys model
                $$r_cur_cm_file = $self->{Pselc_cm_file_path};
            }
            elsif  ($trna_isotype eq "SeCe") {    # use euk selcys model
                $$r_cur_cm_file = $self->{Eselc_cm_file_path};
            }            
        }
    }
}

# find anticodon loop & a-codon from coves or cmsearch output

sub find_anticodon {                

    my $self = shift;
    my ($seq, $ss, $undef_anticodon) = @_;
    my ($antiloop_index, $antiloop, $antiloop_len, $antiloop_end, $ac_index, $anticodon, $verify_ac);

    # Match pattern in secondary structure output, 
    # looking for second stem-loop structure ">>>>...<<<<"
    # that should be the anitocodon stem-loop 

    if ($ss =~ /^([>.]+<[<.]+>[>.]*)>([.]{4,})<+/o) {

        # set to index position of first base in anticodon loop
        $antiloop_index = length($1) + 1;
        $antiloop_len = length($2);   # anticodon loop length
    
        # index of end of anticodon loop
        $antiloop_end = $antiloop_index + $antiloop_len - 1;
    
        $antiloop = substr($seq, $antiloop_index, $antiloop_len);
    
        # remove '-' gaps from loop
        $antiloop =~ s/[\-]//g;      
        # remove introns & non-canonical bases
        $antiloop =~ s/[a-z]//g;      

        # Don't guess if even number of bp in 
        # anticodon loop
        if ((length($antiloop) < 5) || ((length($antiloop) % 2) == 0)) {
            return ($undef_anticodon, -1, -1, -1);
        }
        # get anticodon 
        $ac_index = (length($antiloop) - 3) / 2;
        $anticodon = substr($antiloop, $ac_index, 3);
        $verify_ac = substr($seq, $ac_index + $antiloop_index, 3);
    
        # check to see if anticodon extracted from the entire
        #  trna sequence (coveseq) is same as that extracted from
        #  just the anticodon loop sequence (antiloop)
    
        if ($verify_ac ne $anticodon) {
            return ($undef_anticodon, -1, -1, -1);            
        }
        return ($anticodon, $antiloop_index, $antiloop_end, $ac_index + $antiloop_index + 1);
    }
    else  {
        return ($undef_anticodon, -1, -1, -1);
    }
}

sub find_intron {

    my $self = shift;
    my ($trna_seq, $antiloop_index, $antiloop_end) = @_;
    my ($intron, $istart, $iend, $tmpstr, $antiloop_seq);
    my $min_intron_length = $self->{min_intron_length};

    # check to see if it was unable 
    # to determine the anticodon loop
    if ($antiloop_index == -1) {
        return(0, 0, 0);
    }
    # get subsequence from start of anticodon loop
    # to end of anticodon loop -- look for intron in it
    $antiloop_seq = substr($trna_seq, $antiloop_index, $antiloop_end - $antiloop_index + 1);
    
    if ($antiloop_seq =~ /^(.*[^a-z]+)([a-z]{$min_intron_length,})[^a-z]+/o)  {
        $intron = $2;
    
        # make sure to get the base index for the last (not nec. only) occurrence
        # of the intron sequence string up to end of anticodon loop
        $tmpstr = substr($trna_seq, 0, $antiloop_end+1);
        $istart = index($tmpstr, $intron) + 1; 
        $iend = length($intron) + $istart - 1;
    }
    else {
            $intron = 0; 
    }
    return ($intron, $istart, $iend);
}                        

# is_pseudo_gene
#
# Runs a covariance model without secondary structure 
# information on predicted tRNA, puts this value
# in "hmm_score".  
# Contribution to total score from secondary structure 
# derived by subtracting hmm_score from total score
# Returns non-zero if tRNA scores fall below minima
# for either primary or secondary structure components
# of score

sub is_pseudo_gene {
    
    my $self = shift;
    my $opts = shift;
    my ($r_hmm_score, $r_ss_score, $score, $tmp_trnaseq_file, $trna_name, $tRNA_len) = @_;
    my ($dummy1, $dummy2, $hit_start, $hit_end, $hit_ct, $min_pseudo_filter_score);

    $$r_ss_score = $$r_hmm_score = -1000;   # clear values to be returned
    $dummy1 = $dummy2 = "";                 # return values not used

    # skip check for pseudo gene if score is above minimum
    # -D (disable pseudogene checking) is specified 
    # AND -H option (get hmm scores) is NOT specified
    if ($self->cove_mode()) {
        $min_pseudo_filter_score = $self->{min_cove_pseudo_filter_score};
    }
    elsif ($self->infernal_mode()) {
        $min_pseudo_filter_score = $self->{min_cmsearch_pseudo_filter_score};
    }
    
    if ((($score >= $min_pseudo_filter_score) || $self->{skip_pseudo_filter}) && !$self->{get_hmm_score}) {
        return 0;
    }

    if ($self->cove_mode())
    {
        ($dummy1, $dummy2, $$r_hmm_score) = 
            $self->run_coves($tmp_trnaseq_file, $trna_name, $self->{mainNS_cm_file_path});
    }
    elsif ($self->infernal_mode()) 
    {
       ($$r_hmm_score, $hit_start, $hit_end, $hit_ct) = 
            $self->cmsearch_bestscore($opts, $tmp_trnaseq_file, $trna_name, $tRNA_len, $self->{mainNS_cm_file_path});
    }
    else {
        return -1;                              # Error - no second pass scanner selected
    }
    $$r_ss_score = $score - $$r_hmm_score;      # calc secondary structure
                                                # contribution to total bit score

    if ((($$r_ss_score < $self->{min_ss_score}) || ($$r_hmm_score < $self->{min_hmm_score})) &&
        ($score < $min_pseudo_filter_score)) {
            return 1;
    }
}    

# Get tRNA anticodon, isotype, intron location, and pseudogene status

sub decode_tRNA_properties {

    my $self = shift;
    my $opts = shift;
    my $gc = shift;
    my $log = shift;
    my ($trna_score, $trna_seq, $trna_ss, $r_prescan_tRNA, $trna_start, $trna_end, $cur_cm_file, $tmp_trnaseq_file) = @_;

    my ($anticodon, $acodon_index, $trna_type, $intron, $istart, $iend, @introns,
          $hmm_score, $ss_score, $pseudo_gene_flag,
          $antiloop_index, $antiloop_end, $trna_len, $scan_len);

    $anticodon = "ERR";

    if ($self->cove_mode() || $self->infernal_mode()) {
        ($anticodon, $antiloop_index, $antiloop_end, $acodon_index) = 
            $self->find_anticodon($trna_seq, $trna_ss, $gc->undef_anticodon()); 
    }
    else {
        die "Second pass mode not selected -- can't decode tRNA type\n\n";
    }
    
    # check for problem parsing anticodon loop 
    if (($anticodon eq $gc->undef_anticodon()) || ($trna_seq  eq 'Error'))    
    {
        $anticodon = $gc->undef_anticodon();
        $trna_type = $gc->undef_isotype();
        $intron = 0;        
        
        if ($opts->save_odd_struct()) {     
            open(ODDTRNA, ">>".$opts->odd_struct_file()) ||
                die "FATAL: Can't open ".$opts->odd_struct_file()." to save seconary structures\n\n"; 
            print ODDTRNA "$r_prescan_tRNA->{name} ($trna_start-$trna_end):\n$trna_seq\n$trna_ss\n\n"; 
            close(ODDTRNA);
        }
    }
    else {                               # continue tRNA struct parsing
        ($intron, $istart, $iend) = 
            $self->find_intron($trna_seq, $antiloop_index, $antiloop_end);
        
        if ($intron) {
            push(@introns, {seq=>$intron, start=>$istart, end=>$iend, type=>"CI"});
        }
        
        if (defined $r_prescan_tRNA->{acodon}) {
            if (($anticodon ne (uc($r_prescan_tRNA->{acodon}))) && 
                ($opts->tscan_mode() || $opts->eufind_mode()) && ($opts->strict_params())) {
                $log->write_line("\n$r_prescan_tRNA->{name} - anticondon conflict\t".$opts->second_pass_label().": $anticodon\tfirstpass ($r_prescan_tRNA->{hit_source})".
                    ": $r_prescan_tRNA->{acodon}\n$trna_seq\n$trna_ss\n"); 
            }
        }
        
        $trna_type = $gc->get_tRNA_type($self, $anticodon, $cur_cm_file);
    }
        
    $pseudo_gene_flag = 0;
    $hmm_score = $ss_score = 0;
    
    # Write current tRNA to temp file for re-analysis with other models
    $trna_len = length($trna_seq);
    &write_tRNA($tmp_trnaseq_file, $r_prescan_tRNA->{name}, "", $trna_seq, 1);
    
    if (($trna_type !~ /SeC/) &&
        ($self->is_pseudo_gene($opts, \$hmm_score, \$ss_score, $trna_score, $tmp_trnaseq_file, $r_prescan_tRNA->{name}, $trna_len)) &&
        (!$self->{skip_pseudo_filter}))
    {
        $pseudo_gene_flag = 1;     # set to non-zero for likely
    }                              #  pseudogenes
    
    return ($anticodon, $acodon_index, $trna_type, \@introns, 
            $hmm_score, $ss_score, $pseudo_gene_flag);
}


sub scan_split_tRNAs {

    my $self = shift;
    my $opts = shift;
    my $constants = shift;
    my $stats = shift;
    my $gc = shift;
    my $log = shift;
    my $r_sec_pass_hits = shift;
    
    my $r_pair = {};
    my ($r_pairs_index, $r_five_half_hits, $r_three_half_hits) = $self->scan_split_tRNA_halves($opts, $constants, $stats, $gc, $log, $r_sec_pass_hits);   
    my ($r_pairs) = $self->scan_split_tRNAs_in_long_introns($opts, $constants, $stats, $gc, $log, $r_sec_pass_hits);
    
    if (scalar(@$r_pairs) > 0) {
        my $five_half_count = scalar(@$r_five_half_hits);
        my $three_half_count = scalar(@$r_three_half_hits);
        
        for (my $ct = 0; $ct < scalar(@$r_pairs); $ct++) {
            push(@$r_five_half_hits, $r_pairs->[$ct]->{"5h"});
            push(@$r_three_half_hits, $r_pairs->[$ct]->{"3h"});
            $r_pair = {};
            $r_pair->{"5h"} = $five_half_count;
            $r_pair->{"3h"} = $three_half_count;
            push(@$r_pairs_index, $r_pair);
            $five_half_count++;
            $three_half_count++;
        }
    }
    
    &output_split_fragments($opts, $r_pairs_index, $r_five_half_hits, $r_three_half_hits);
}

sub scan_split_tRNA_halves {

    my $self = shift;
    my $opts = shift;
    my $constants = shift;
    my $stats = shift;
    my $gc = shift;
    my $log = shift;
    my $r_sec_pass_hits = shift;
    
    my $tmp_trnaseq_file = $constants->tmp_trnaseq_file();
    
    my (@five_half_hits, @three_half_hits, $five_half_count, $three_half_count, $r_five_half_hit, $r_three_half_hit, $index_5h, $index_3h,
        $acceptor_half_seq, $rc_seq, $rc_acceptor_half_seq, $r_pair, @pairs);
    
    my ($r_valid, $match, $intron_len, $r_intron, $split_trna_ct);
    
    my ($cur_cm_file, $cms_output, @half_hit_list, $r_cm_hit);
    
    $split_trna_ct = 0;
    
    my $seq_file = tRNAscanSE::Sequence->new;
    my @sorted_cm_hits = sort sort_cm_hits_by_start @$r_sec_pass_hits;
    
    # find 5'half and 3'half tRNA hits using cmsearch
    $seq_file->mask_out_sequence($opts->fasta_file(), $constants->tmp_masked_fa(), \@sorted_cm_hits);
    $cur_cm_file = $self->{arch_five_half_cm_file_path};
    if (!$self->run_gw_cmsearch($opts, $log, \@five_half_hits, \$cur_cm_file, $constants->tmp_masked_fa(), "", 1)) {
        return 0;
    }
    $cur_cm_file = $self->{arch_three_half_cm_file_path};
    if (!$self->run_gw_cmsearch($opts, $log, \@three_half_hits, \$cur_cm_file, $constants->tmp_masked_fa(), "", 1)) {
        return 0;
    }

    foreach $r_cm_hit (@$r_sec_pass_hits) {
        $intron_len = 0;
        $match = 0;
        if (scalar(@{$r_cm_hit->{introns}}) > 0) {           
            foreach $r_intron (@{$r_cm_hit->{introns}}) {
                if ($r_intron->{type} eq "CI") {
                    $intron_len = length($r_intron->{seq});
                }
            }
        }
        if ($r_cm_hit->{score} < $self->{split_tRNA_scan_cutoff}) {
            $r_valid = &valid_structure($r_cm_hit->{ss}, $intron_len);
            if (!$r_valid->{tRNA}) {
                @half_hit_list = ();
                &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $r_cm_hit->{seq}, 1);
                $cur_cm_file = $self->{arch_five_half_cm_file_path};
                if (!$self->exec_cmsearch(0, \$cur_cm_file, 0, $tmp_trnaseq_file, $r_cm_hit->{seqname}, \$cms_output)) {
                    return 0;
                }
                $self->parse_cmsearch($cms_output, \@half_hit_list, 0, $r_cm_hit, 1);
                foreach $r_five_half_hit (@half_hit_list) {
                    if ($r_five_half_hit->{score} >= $self->{half_tRNA_cutoff}) {
                        push(@five_half_hits, $r_five_half_hit);
                        $match = 1;
                        last;
                    }
                }
                
                if (!$match) {
                    @half_hit_list = ();
                    $cur_cm_file = $self->{arch_three_half_cm_file_path};
                    if (!$self->exec_cmsearch(0, \$cur_cm_file, 0, $tmp_trnaseq_file, $r_cm_hit->{seqname}, \$cms_output)) {
                        return 0;
                    }
                    $self->parse_cmsearch($cms_output, \@half_hit_list, 0, $r_cm_hit, 1);
                    foreach $r_three_half_hit (@half_hit_list) {
                        if ($r_three_half_hit->{score} >= $self->{half_tRNA_cutoff}) {
                            push(@three_half_hits, $r_three_half_hit);
                            $match = 1;
                            last;
                        }
                    }                    
                }
            }
        }
    }

    $five_half_count = 0; $three_half_count = 0;
    foreach $r_five_half_hit (@five_half_hits) {
        $r_five_half_hit->{seqname} = $r_five_half_hit->{hit_seqname};
        $r_five_half_hit->{start} = $r_five_half_hit->{tRNA_start};
        $r_five_half_hit->{end} = $r_five_half_hit->{tRNA_end};
        $five_half_count++ if ($r_five_half_hit->{score} >= $self->{half_tRNA_cutoff});
    }
    foreach $r_three_half_hit (@three_half_hits) {
        $r_three_half_hit->{seqname} = $r_three_half_hit->{hit_seqname};
        $r_three_half_hit->{start} = $r_three_half_hit->{tRNA_start};
        $r_three_half_hit->{end} = $r_three_half_hit->{tRNA_end};
        $three_half_count++ if ($r_three_half_hit->{score} >= $self->{half_tRNA_cutoff});
    }

    @five_half_hits = sort sort_cm_hits_by_start @five_half_hits;
    @three_half_hits = sort sort_cm_hits_by_start @three_half_hits;
    
    if ($five_half_count >= $three_half_count) {
        for ($index_5h = 0; $index_5h < scalar(@five_half_hits); $index_5h++) {
            if ($five_half_hits[$index_5h]->{score} >= $self->{half_tRNA_cutoff}) {
                $acceptor_half_seq = uc(&get_acceptor_half($five_half_hits[$index_5h]->{seq}, $five_half_hits[$index_5h]->{ss}, "5h"));
                $rc_seq = &rev_comp_seq($acceptor_half_seq);
                $rc_seq =~ s/C/\[CT\]/g;
                foreach ($index_3h = 0; $index_3h < scalar(@three_half_hits); $index_3h++) {
                    if ($three_half_hits[$index_3h]->{score} >= $self->{half_tRNA_cutoff}) {
                        $rc_acceptor_half_seq = uc(&get_acceptor_half($three_half_hits[$index_3h]->{seq}, $three_half_hits[$index_3h]->{ss}, "3h"));
                        if ($rc_acceptor_half_seq =~ /$rc_seq/) {
                            $r_pair = {};
                            $r_pair->{"5h"} = $index_5h;
                            $r_pair->{"3h"} = $index_3h;
                            push(@pairs, $r_pair);
                            $five_half_hits[$index_5h]->{pair} = 1;
                            $three_half_hits[$index_3h]->{pair} = 1;
                        }
                    }
                }
            }
        }
    }
    else {
        for ($index_3h = 0; $index_3h < scalar(@three_half_hits); $index_3h++) {
            if ($three_half_hits[$index_3h]->{score} >= $self->{half_tRNA_cutoff}) {
                $acceptor_half_seq = uc(&get_acceptor_half($three_half_hits[$index_3h]->{seq}, $three_half_hits[$index_3h]->{ss}, "3h"));
                $rc_seq = &rev_comp_seq($acceptor_half_seq);
                $rc_seq =~ s/C/\[CT\]/g;
                foreach ($index_5h = 0; $index_5h < scalar(@five_half_hits); $index_5h++) {
                    if ($five_half_hits[$index_5h]->{score} >= $self->{half_tRNA_cutoff}) {
                        $rc_acceptor_half_seq = uc(&get_acceptor_half($five_half_hits[$index_5h]->{seq}, $five_half_hits[$index_5h]->{ss}, "5h"));
                        if ($rc_acceptor_half_seq =~ /$rc_seq/) {
                            $r_pair = {};
                            $r_pair->{"5h"} = $index_5h;
                            $r_pair->{"3h"} = $index_3h;
                            push(@pairs, $r_pair);
                            $five_half_hits[$index_5h]->{pair} = 1;
                            $three_half_hits[$index_3h]->{pair} = 1;
                        }
                    }
                }
            }
        }        
    }

    for ($index_5h = 0; $index_5h < scalar(@five_half_hits); $index_5h++) {
        if ($five_half_hits[$index_5h]->{score} >= $self->{half_tRNA_cutoff}) {
            if (!defined $five_half_hits[$index_5h]->{pair}) {
                $r_pair = {};
                $r_pair->{"5h"} = $index_5h;
                push(@pairs, $r_pair);
            }
        }
    }
    for ($index_3h = 0; $index_3h < scalar(@three_half_hits); $index_3h++) {
        if ($three_half_hits[$index_3h]->{score} >= $self->{half_tRNA_cutoff}) {
            if (!defined $three_half_hits[$index_3h]->{pair}) {
                $r_pair = {};
                $r_pair->{"3h"} = $index_3h;
                push(@pairs, $r_pair);
            }
        }
    }
    
    return (\@pairs, \@five_half_hits, \@three_half_hits);
}

sub scan_split_tRNAs_in_long_introns {

    my $self = shift;
    my $opts = shift;
    my $constants = shift;
    my $stats = shift;
    my $gc = shift;
    my $log = shift;
    my $r_sec_pass_hits = shift;
    
    my $tmp_trnaseq_file = $constants->tmp_trnaseq_file();
    my ($r_valid, $scan_flag, $intron_len, $skip, $scan_trna_seq, $over_cutoff_count,
        @rescan_trna_hits, @five_half_hit_list, @three_half_hit_list, @pairs, $r_pair, $trna_ct, $split_trna_ct);
    
    my ($cur_cm_file, $cms_output, $r_cm_hit, $r_cm_hit2, $r_cm_hit_temp, $r_intron, @intron_hit_list, $intron_idx);
    
    $split_trna_ct = 0;
    foreach $r_cm_hit (@$r_sec_pass_hits) {
        $scan_flag = 0;
        $intron_len = 0;
        $intron_idx = -1;
        if (scalar(@{$r_cm_hit->{introns}}) > 0) {           
            foreach $r_intron (@{$r_cm_hit->{introns}}) {
                $intron_idx++;
                if ((length($r_intron->{seq}) > 100) && ($r_cm_hit->{isotype} ne "Trp") && ($r_cm_hit->{isotype} ne "Tyr")) {
                    $scan_flag = 1;
                    last;
                }
            }
        }
        if ($scan_flag) {
            $r_cm_hit_temp = {};
            $r_cm_hit_temp->{start} = 1;
            $r_cm_hit_temp->{strand} = 1;
            $r_cm_hit_temp->{hit_source} = $r_cm_hit->{hit_source};
            $r_cm_hit_temp->{seqname} = $r_cm_hit->{seqname};
            $r_cm_hit_temp->{src_seqname} = $r_cm_hit->{seqname};
            $r_cm_hit_temp->{upstream} = "";
            $r_cm_hit_temp->{downstream} = "";

             $scan_trna_seq = uc(substr($r_cm_hit->{seq}, $r_cm_hit->{introns}->[$intron_idx]->{start} - 13, $self->{left_splicing_len}) .
                substr($r_cm_hit->{seq}, $r_cm_hit->{introns}->[$intron_idx]->{end} - $self->{right_splicing_len} + 8, $self->{right_splicing_len}));
            
            &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $scan_trna_seq, 1);
            
            $cur_cm_file = $self->{arch_intron_cm_file_path};
            if (!$self->run_cmsearch_intron($opts, $log, \@intron_hit_list, \$cur_cm_file, $r_cm_hit, $tmp_trnaseq_file, \$over_cutoff_count)) {
                return 0;
            }
            if ($over_cutoff_count > 0) {
                $scan_trna_seq = $r_cm_hit->{seq};
                $scan_trna_seq =~ s/$r_cm_hit->{introns}->[$intron_idx]->{seq}//;
                $r_cm_hit_temp->{seq} = $scan_trna_seq;
                $r_cm_hit_temp->{len} = length($scan_trna_seq);
                $r_cm_hit_temp->{end} = $r_cm_hit_temp->{len};
                &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $scan_trna_seq, 1);

                if ($opts->cove_mode()) 
                {
                    $trna_ct = $self->analyze_with_cove($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                               $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                }
                elsif ($opts->infernal_mode()) 
                {
                    $trna_ct = $self->analyze_with_cmsearch($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                                   $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                }
                
                if ((scalar(@rescan_trna_hits) > 0) && ($rescan_trna_hits[0]->{score} >= $r_cm_hit->{score})) {
                    $split_trna_ct++;

                    &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $scan_trna_seq, 1);

                    $cur_cm_file = $self->{arch_five_half_cm_file_path};
                    if (!$self->exec_cmsearch(0, \$cur_cm_file, 0, $tmp_trnaseq_file, $r_cm_hit->{seqname}, \$cms_output)) {
                        return 0;
                    }
                    $self->parse_cmsearch($cms_output, \@five_half_hit_list, 0, $rescan_trna_hits[0], 1);

                    $cur_cm_file = $self->{arch_three_half_cm_file_path};
                    if (!$self->exec_cmsearch(0, \$cur_cm_file, 0, $tmp_trnaseq_file, $r_cm_hit->{seqname}, \$cms_output)) {
                        return 0;
                    }
                    $self->parse_cmsearch($cms_output, \@three_half_hit_list, 0, $rescan_trna_hits[0], 1);

                    if (scalar(@five_half_hit_list) > 0 && scalar(@three_half_hit_list) > 0) {
                        $five_half_hit_list[0]->{tRNA_start} = $r_cm_hit->{start};
                        if ($r_cm_hit->{strand}) {	
                            $five_half_hit_list[0]->{tRNA_end} = ($r_cm_hit->{introns}->[$intron_idx]->{start} + $r_cm_hit->{start} - 2);
    						$three_half_hit_list[0]->{tRNA_start} = ($r_cm_hit->{introns}->[$intron_idx]->{end} + $r_cm_hit->{start});
                        }
                        else {
                            $five_half_hit_list[0]->{tRNA_end} = ($r_cm_hit->{start} - $r_cm_hit->{introns}->[$intron_idx]->{start} + 2);
                            $three_half_hit_list[0]->{tRNA_start} = ($r_cm_hit->{start} - $r_cm_hit->{introns}->[$intron_idx]->{end});
                        }                        
                        $three_half_hit_list[0]->{tRNA_end} = $r_cm_hit->{end};
                        $r_pair = {};
                        $r_pair->{"5h"} = $five_half_hit_list[0];
                        $r_pair->{"3h"} = $three_half_hit_list[0];
                        push(@pairs, $r_pair);
                    }
                }
            }
        }
    }
    return \@pairs;
}

sub scan_noncanonical_introns {

    my $self = shift;
    my $opts = shift;
    my $constants = shift;
    my $stats = shift;
    my $gc = shift;
    my $log = shift;
    my $seq_file = shift;
    my $r_sec_pass_hits = shift;
    
    my $tmp_trnaseq_file = $constants->tmp_trnaseq_file();
    my ($r_valid, $scan_flag, $skip, $ci_intron_index, $scan_trna_seq, $over_cutoff_count,
        $best_score, $total_intron_len, $new_start, $r_new_intron, @rescan_trna_hits, $trna_ct);
    my ($partial_scan_trna_seq, @partial_intron_hit_list, $partial_over_cutoff_count, $last_end);
    
    my ($cur_cm_file, $r_cm_hit, $r_cm_hit_temp, @extra_cm_hit_list, $r_intron, $r_intron_hit, @intron_hit_list, $tRNA_seq, $upstream, $downstream);
    my ($anticodon, $acodon_index, $isotype, $r_introns, $hmm_score, $ss_score, $pseudo_gene_flag);
    my ($cur_tRNA_ct) = 0;
    
    my $masked_seq_file = tRNAscanSE::Sequence->new;
    my @sorted_cm_hits = sort sort_cm_hits_by_start @$r_sec_pass_hits;
   
    # find extra tRNA hits using cmsearch
    $masked_seq_file->mask_out_sequence($opts->fasta_file(), $constants->tmp_masked_fa(), \@sorted_cm_hits);
    if (!$self->run_gw_cmsearch($opts, $log, \@extra_cm_hit_list, \$cur_cm_file, $constants->tmp_masked_fa(), "", 0)) {
        return 0;
    }

    foreach my $r_extra_hit (@extra_cm_hit_list)
    {
        if ($r_extra_hit->{score} >= $self->{cm_cutoff}) {
            ($tRNA_seq, $upstream, $downstream) = $seq_file->get_tRNA_sequence($r_extra_hit->{hit_seqname}, $r_extra_hit->{strand},
                                                                               $r_extra_hit->{tRNA_start}, $r_extra_hit->{tRNA_end},
                                                                               $log, $opts, $constants);
            $r_extra_hit->{name} = $r_extra_hit->{hit_seqname};
            
            &write_tRNA($tmp_trnaseq_file, $r_extra_hit->{hit_seqname}, " ", $r_extra_hit->{seq}, 1);
    
            ($anticodon, $acodon_index, $isotype, $r_introns, $hmm_score, $ss_score, $pseudo_gene_flag) = 
                 $self->decode_tRNA_properties ($opts, $gc, $log, $r_extra_hit->{score}, $r_extra_hit->{seq}, $r_extra_hit->{ss}, $r_extra_hit,
                              $r_extra_hit->{tRNA_start}, $r_extra_hit->{tRNA_end}, $cur_cm_file, $tmp_trnaseq_file);
            
            $r_cm_hit = {};     
            $r_cm_hit =
                {seqname =>$r_extra_hit->{hit_seqname}, score=>$r_extra_hit->{score}, ss=>$r_extra_hit->{ss}, seq=>$r_extra_hit->{seq}, model=>$r_extra_hit->{model},
                 start=>$r_extra_hit->{tRNA_start}, end=>$r_extra_hit->{tRNA_end}, len=>$r_extra_hit->{tRNA_len}, ID=>$r_extra_hit->{name},
                 acodon=>$anticodon, acodon_pos =>$acodon_index, isotype=>$isotype,
                 introns=>$r_introns, hmm_score=>$hmm_score, 
                 ss_score=>$ss_score, is_pseudo=>$pseudo_gene_flag,
                 src_seqlen=>3, src_seqname=>$r_extra_hit->{hit_seqname},
                 strand=>$r_extra_hit->{strand}, hit_source=>"Inf",
                 upstream=>$upstream, downstream=>$downstream, extra=>1};
                
            push (@$r_sec_pass_hits, $r_cm_hit);
        }
    }
    
    # scan for noncanonical introns
    $cur_cm_file = $self->{arch_intron_cm_file_path};
    foreach $r_cm_hit (sort sort_cm_hits_for_output @$r_sec_pass_hits) {
        
        $scan_flag = 0;
        
        # renumber tRNA hits
        $cur_tRNA_ct++;
        $r_cm_hit->{ID} = $r_cm_hit->{seqname}.".t".$cur_tRNA_ct;
        
        if ($r_cm_hit->{extra}) {
            $scan_flag = 1;
        }
        else {
            if (scalar(@{$r_cm_hit->{introns}}) > 0) {
                $r_valid = &valid_structure($r_cm_hit->{ss}, length($r_cm_hit->{introns}->[0]->{seq}));
            }
            else {
                $r_valid = &valid_structure($r_cm_hit->{ss}, 0);
            }
            if (scalar(@{$r_cm_hit->{introns}}) > 0) {
                $scan_flag = 1;
            }
            elsif ($r_cm_hit->{acodon} eq "???") {
                $scan_flag = 1;
            }
            elsif ($r_cm_hit->{score} < $self->{nci_scan_cutoff}) {
                if (!$r_valid->{tRNA}) {
                    $scan_flag = 1;
                }
            }
        }
        
        if ($scan_flag) {

            $log->write_line("Scan for noncanonical intron ".$r_cm_hit->{ID}."-".$r_cm_hit->{isotype}.$r_cm_hit->{acodon});
 
            $total_intron_len = 0;
            @intron_hit_list = ();
            $best_score = $r_cm_hit->{score};
            $r_cm_hit_temp = {};
            $r_cm_hit_temp->{start} = 1;
            $r_cm_hit_temp->{strand} = 1;
            $r_cm_hit_temp->{hit_source} = $r_cm_hit->{hit_source};
            $r_cm_hit_temp->{seqname} = $r_cm_hit->{seqname};
            $r_cm_hit_temp->{src_seqname} = $r_cm_hit->{seqname};
            $r_cm_hit_temp->{upstream} = "";
            $r_cm_hit_temp->{downstream} = "";
            
            $scan_trna_seq = uc($r_cm_hit->{upstream}.$r_cm_hit->{seq}.$r_cm_hit->{downstream});
            $partial_scan_trna_seq = uc($r_cm_hit->{upstream}.substr($r_cm_hit->{seq}, 0, 12));
            $r_cm_hit_temp->{seq} = $scan_trna_seq;
            $r_cm_hit_temp->{precursor} = $scan_trna_seq;
            $r_cm_hit_temp->{len} = length($scan_trna_seq);
            $r_cm_hit_temp->{end} = $r_cm_hit_temp->{len};
            
            &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $scan_trna_seq, 1);
            
            if (!$self->run_cmsearch_intron($opts, $log, \@intron_hit_list, \$cur_cm_file, $r_cm_hit, $tmp_trnaseq_file, \$over_cutoff_count)) {
                return 0;
            }

            &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $partial_scan_trna_seq, 1);
            
            if (!$self->run_cmsearch_intron($opts, $log, \@partial_intron_hit_list, \$cur_cm_file, $r_cm_hit, $tmp_trnaseq_file, \$partial_over_cutoff_count)) {
                return 0;
            }
            foreach $r_intron_hit (@intron_hit_list) {
                $skip = 0;
                $ci_intron_index = -1;
                $r_intron_hit->{overlap} = $ci_intron_index;
                if ($r_intron_hit->{score} >= $self->{BHB_cm_cutoff}) {
                    foreach $r_intron (@{$r_cm_hit->{introns}}) {
                        $ci_intron_index++;
                        if ($r_intron->{type} eq "CI") {
                            if ($r_intron_hit->{intron_seq} eq uc($r_intron->{seq})) {
                                $skip = 1;
                                last;
                            }
                            elsif (($r_intron_hit->{start} >= $r_intron->{start} && $r_intron_hit->{start} <= $r_intron->{end}) ||
                                    ($r_intron_hit->{end} >= $r_intron->{start} && $r_intron_hit->{end} <= $r_intron->{end})) {
                                $r_intron_hit->{overlap} = $ci_intron_index;
                                last;
                            }
                        }
                    }
                    if (!$skip) {
                        if (($r_intron_hit->{tRNA_start} < $r_cm_hit->{start} && $r_cm_hit->{strand}) ||
                            ($r_intron_hit->{tRNA_start} > $r_cm_hit->{start} && !$r_cm_hit->{strand})) {
                            if ($r_valid->{acceptor}) {
                                $skip = 1;
                            }
                        }
                    }
                    if (!$skip) {
                        @rescan_trna_hits = ();
                        my $idx_intron = index($scan_trna_seq, $r_intron_hit->{intron_seq});
                        $scan_trna_seq = substr($scan_trna_seq, 0, $idx_intron) . substr($scan_trna_seq, $idx_intron + $r_intron_hit->{intron_len});
                        $r_cm_hit_temp->{seq} = $scan_trna_seq;
                        $r_cm_hit_temp->{len} = length($scan_trna_seq);
                        $r_cm_hit_temp->{end} = $r_cm_hit_temp->{len};
                        &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $scan_trna_seq, 1);
                        
                        if ($opts->cove_mode()) 
                        {
                            $trna_ct = $self->analyze_with_cove($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                                       $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                        }
                        elsif ($opts->infernal_mode()) 
                        {
                            $trna_ct = $self->analyze_with_cmsearch($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                                           $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                        }
                        
                        if (scalar(@rescan_trna_hits) > 0) {
                            if ($rescan_trna_hits[0]->{score} > $best_score) {
                                $best_score = $rescan_trna_hits[0]->{score};
                                $total_intron_len += $r_intron_hit->{intron_len};
                                $r_cm_hit->{acodon} = $rescan_trna_hits[0]->{acodon};
                                $r_cm_hit->{acodon_pos} = $rescan_trna_hits[0]->{acodon_pos};
                                $r_cm_hit->{isotype} = $rescan_trna_hits[0]->{isotype};
                                $r_cm_hit->{score} = $rescan_trna_hits[0]->{score};
                                $r_cm_hit->{hmm_score} = $rescan_trna_hits[0]->{hmm_score};
                                $r_cm_hit->{ss_score} = $rescan_trna_hits[0]->{ss_score};
                                $r_cm_hit->{is_pseudo} = $rescan_trna_hits[0]->{is_pseudo};
                                $r_cm_hit->{ss} = $rescan_trna_hits[0]->{ss};
                                $r_cm_hit->{seq} = $rescan_trna_hits[0]->{seq};
                                $r_cm_hit->{precursor} = substr($r_cm_hit_temp->{precursor}, $rescan_trna_hits[0]->{start} - 1, $rescan_trna_hits[0]->{len} + $total_intron_len);
                                $r_cm_hit->{model} = $rescan_trna_hits[0]->{model};
                                $r_cm_hit->{len} = $rescan_trna_hits[0]->{len} + $total_intron_len;
                                $new_start = $rescan_trna_hits[0]->{start};
                                if ($r_intron_hit->{overlap} > -1) {
                                    $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{seq} = $r_intron_hit->{intron_seq};
                                    $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{start} = $r_intron_hit->{start};
                                    $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{end} = $r_intron_hit->{end};
                                    $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{type} = "NCI";
                                }
                                else {
                                    push(@{$r_cm_hit->{introns}}, {seq=>$r_intron_hit->{intron_seq}, start=>$r_intron_hit->{start},
                                                                   end=>$r_intron_hit->{end}, type=>"NCI"});
                                }
                                if (scalar(@{$rescan_trna_hits[0]->{introns}}) > 0) {
                                    $r_new_intron = $rescan_trna_hits[0]->{introns}->[0];
                                    my $set_ci = 0;
                                    foreach $r_intron (@{$r_cm_hit->{introns}}) {
                                        if ($r_intron->{type} eq "CI") {
                                            $r_intron = $r_new_intron;
                                            $set_ci = 1;
                                        }
                                    }
                                    if (!$set_ci) {
                                        push(@{$r_cm_hit->{introns}}, {seq=>$r_new_intron->{seq}, start=>$r_new_intron->{start},
                                                                       end=>$r_new_intron->{end}, type=>"CI"});                                        
                                    }
                                }
                                else {
                                    foreach $r_intron (@{$r_cm_hit->{introns}}) {
                                        if ($r_intron->{type} eq "CI") {
                                            $r_intron = undef;
                                        }
                                    }                                    
                                }
                            }
                        }
                    }
                }
            }
            if ($total_intron_len > 0) {
                if ($r_cm_hit->{strand}) {
                    $r_cm_hit->{start} = $r_cm_hit->{start} - length($r_cm_hit->{upstream}) + $new_start - 1;
                    $r_cm_hit->{end} = $r_cm_hit->{start} + $r_cm_hit->{len} - 1;
                }
                else {
                    $r_cm_hit->{start} = $r_cm_hit->{start} + length($r_cm_hit->{upstream}) - $new_start + 1;
                    $r_cm_hit->{end} = $r_cm_hit->{start} - $r_cm_hit->{len} + 1;                                        
                }
                foreach $r_intron (@{$r_cm_hit->{introns}}) {
                    if (defined $r_intron) {
                        $r_intron->{start} = index($r_cm_hit->{precursor}, uc($r_intron->{seq})) + 1;
                        $r_intron->{end} = $r_intron->{start} + length($r_intron->{seq}) - 1;
                        if (($r_intron->{start} == ($r_cm_hit->{acodon_pos}+4)) && ($r_intron->{type} eq "NCI")) {
                            $r_intron->{type} = "CI";
                        }
                    }
                }
                @{$r_cm_hit->{introns}} = sort sort_intron_by_start @{$r_cm_hit->{introns}};
                $last_end = -1;
                for (my $ct = 0; $ct < scalar(@{$r_cm_hit->{introns}}); $ct++) {
                    if ($last_end == $r_cm_hit->{introns}->[$ct]->{start} - 1) {
                        $last_end = $r_cm_hit->{introns}->[$ct]->{end};
                        $r_cm_hit->{introns}->[$ct]->{start} = $r_cm_hit->{introns}->[$ct - 1]->{start};
                        $r_cm_hit->{introns}->[$ct]->{seq} = $r_cm_hit->{introns}->[$ct - 1]->{seq} . $r_cm_hit->{introns}->[$ct]->{seq};
                        $r_cm_hit->{introns}->[$ct - 1] = undef;
                    }
                    $last_end = $r_cm_hit->{introns}->[$ct]->{end};
                }
            }
            elsif(!$r_valid->{acceptor} && $r_cm_hit->{score} > $self->{split_tRNA_scan_cutoff}) {
                foreach $r_intron_hit (@partial_intron_hit_list) {
                    $skip = 0;
                    $ci_intron_index = -1;
                    $r_intron_hit->{overlap} = $ci_intron_index;
                    if ($r_intron_hit->{score} >= $self->{BHB_cm_cutoff}) {
                        @rescan_trna_hits = ();
                        my $idx_intron = index($partial_scan_trna_seq, $r_intron_hit->{intron_seq});
                        $partial_scan_trna_seq = substr($partial_scan_trna_seq, 0, $idx_intron) . substr($partial_scan_trna_seq, $idx_intron + $r_intron_hit->{intron_len}) .
                            substr($r_cm_hit->{seq}, 12) . $r_cm_hit->{downstream};
                        $r_cm_hit_temp->{seq} = $partial_scan_trna_seq;
                        $r_cm_hit_temp->{len} = length($partial_scan_trna_seq);
                        $r_cm_hit_temp->{end} = $r_cm_hit_temp->{len};
                        &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $partial_scan_trna_seq, 1);
                        
                        if ($opts->cove_mode()) 
                        {
                            $trna_ct = $self->analyze_with_cove($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                                       $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                        }
                        elsif ($opts->infernal_mode()) 
                        {
                            $trna_ct = $self->analyze_with_cmsearch($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                                           $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                        }
                        
                        if (scalar(@rescan_trna_hits) > 0) {
                            if ($rescan_trna_hits[0]->{score} > $best_score) {
                                $best_score = $rescan_trna_hits[0]->{score};
                                $total_intron_len += $r_intron_hit->{intron_len};
                                $r_cm_hit->{acodon} = $rescan_trna_hits[0]->{acodon};
                                $r_cm_hit->{acodon_pos} = $rescan_trna_hits[0]->{acodon_pos};
                                $r_cm_hit->{isotype} = $rescan_trna_hits[0]->{isotype};
                                $r_cm_hit->{score} = $rescan_trna_hits[0]->{score};
                                $r_cm_hit->{hmm_score} = $rescan_trna_hits[0]->{hmm_score};
                                $r_cm_hit->{ss_score} = $rescan_trna_hits[0]->{ss_score};
                                $r_cm_hit->{is_pseudo} = $rescan_trna_hits[0]->{is_pseudo};
                                $r_cm_hit->{ss} = $rescan_trna_hits[0]->{ss};
                                $r_cm_hit->{seq} = $rescan_trna_hits[0]->{seq};
                                $r_cm_hit->{precursor} = substr($r_cm_hit_temp->{precursor}, $rescan_trna_hits[0]->{start} - 1, $rescan_trna_hits[0]->{len} + $total_intron_len);
                                $r_cm_hit->{model} = $rescan_trna_hits[0]->{model};
                                $r_cm_hit->{len} = $rescan_trna_hits[0]->{len} + $total_intron_len;
                                $new_start = $rescan_trna_hits[0]->{start};
                                if (scalar(@{$rescan_trna_hits[0]->{introns}}) > 0) {
                                    $r_new_intron = $rescan_trna_hits[0]->{introns}->[0];
                                }
                                push(@{$r_cm_hit->{introns}}, {seq=>$r_intron_hit->{intron_seq}, start=>$r_intron_hit->{start},
                                                               end=>$r_intron_hit->{end}, type=>"NCI"});
                            }
                        }
                    }
                }
                if ($total_intron_len > 0) {
                    if ($r_cm_hit->{strand}) {
                        $r_cm_hit->{start} = $r_cm_hit->{start} - length($r_cm_hit->{upstream}) + $new_start - 1;
                        $r_cm_hit->{end} = $r_cm_hit->{start} + $r_cm_hit->{len} - 1;
                    }
                    else {
                        $r_cm_hit->{start} = $r_cm_hit->{start} + length($r_cm_hit->{upstream}) - $new_start + 1;
                        $r_cm_hit->{end} = $r_cm_hit->{start} - $r_cm_hit->{len} + 1;                                        
                    }
                    foreach $r_intron (@{$r_cm_hit->{introns}}) {
                        if (defined $r_intron) {
    #                        if ((($r_intron->{type} eq "CI") && ($r_intron->{seq} eq $r_new_intron->{seq})) || ($r_intron->{type} eq "NCI")) {
                                $r_intron->{start} = index($r_cm_hit->{precursor}, uc($r_intron->{seq})) + 1;
                                $r_intron->{end} = $r_intron->{start} + length($r_intron->{seq}) - 1;
    #                        }
                            if (($r_intron->{start} == ($r_cm_hit->{acodon_pos}+4)) && ($r_intron->{type} eq "NCI")) {
                                $r_intron->{type} = "CI";
                            }
                        }
                    }
                    @{$r_cm_hit->{introns}} = sort sort_intron_by_start @{$r_cm_hit->{introns}};
                    $last_end = -1;
                    for (my $ct = 0; $ct < scalar(@{$r_cm_hit->{introns}}); $ct++) {
                        if ($last_end == $r_cm_hit->{introns}->[$ct]->{start} - 1) {
                            $last_end = $r_cm_hit->{introns}->[$ct]->{end};
                            $r_cm_hit->{introns}->[$ct]->{start} = $r_cm_hit->{introns}->[$ct - 1]->{start};
                            $r_cm_hit->{introns}->[$ct]->{seq} = $r_cm_hit->{introns}->[$ct - 1]->{seq} . $r_cm_hit->{introns}->[$ct]->{seq};
                            $r_cm_hit->{introns}->[$ct - 1] = undef;
                        }
                        $last_end = $r_cm_hit->{introns}->[$ct]->{end};
                    }
                }                
            }
             if (scalar(@{$r_cm_hit->{introns}}) > 0) {
                if ($r_cm_hit->{introns}->[0]->{type} eq "CI") {
                    $r_valid = &valid_structure($r_cm_hit->{ss}, length($r_cm_hit->{introns}->[0]->{seq}));
                }
                else {
                    $r_valid = &valid_structure($r_cm_hit->{ss}, 0);
                }
            }
            elsif ($r_cm_hit->{score} > $self->{split_tRNA_scan_cutoff}) {
               $r_valid = &valid_structure($r_cm_hit->{ss}, 0);
            }
            if (!$r_valid->{tRNA} && $r_cm_hit->{score} > $self->{split_tRNA_scan_cutoff}) {
                @intron_hit_list = ();
                &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $scan_trna_seq, 1);
                if (!$self->run_cmsearch_intron($opts, $log, \@intron_hit_list, \$cur_cm_file, $r_cm_hit, $tmp_trnaseq_file, \$over_cutoff_count)) {
                    return 0;
                }
                foreach $r_intron_hit (@intron_hit_list) {
                    $skip = 0;
                    $ci_intron_index = -1;
                    $r_intron_hit->{overlap} = $ci_intron_index;
                    if ($r_intron_hit->{score} >= $self->{BHB_cm_cutoff}) {
                        foreach $r_intron (@{$r_cm_hit->{introns}}) {
                            $ci_intron_index++;
                            if ($r_intron->{type} eq "CI") {
                                if ($r_intron_hit->{intron_seq} eq uc($r_intron->{seq})) {
                                    $skip = 1;
                                    last;
                                }
                                else {
                                    my $intron_hit_start = index($r_cm_hit->{precursor}, uc($r_intron_hit->{intron_seq})) + 1;
                                    my $intron_hit_end = $intron_hit_start + length($r_intron_hit->{intron_seq}) - 1;
                                    if (($intron_hit_start >= $r_intron->{start} && $intron_hit_start <= $r_intron->{end}) ||
                                        ($intron_hit_end >= $r_intron->{start} && $intron_hit_end <= $r_intron->{end})) {
                                        $r_intron_hit->{overlap} = $ci_intron_index;
                                        last;
                                    }
                                }
                            }
                        }
                        if (!$skip) {
                            @rescan_trna_hits = ();
                            my $idx_intron = index($scan_trna_seq, $r_intron_hit->{intron_seq});
                            $scan_trna_seq = substr($scan_trna_seq, 0, $idx_intron) . substr($scan_trna_seq, $idx_intron + $r_intron_hit->{intron_len});
                            $r_cm_hit_temp->{seq} = $scan_trna_seq;
                            $r_cm_hit_temp->{len} = length($scan_trna_seq);
                            $r_cm_hit_temp->{end} = $r_cm_hit_temp->{len};
                            &write_tRNA($tmp_trnaseq_file, $r_cm_hit->{seqname}, " ", $scan_trna_seq, 1);
                            
                            if ($opts->cove_mode()) 
                            {
                                $trna_ct = $self->analyze_with_cove($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                                           $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                            }
                            elsif ($opts->infernal_mode()) 
                            {
                                $trna_ct = $self->analyze_with_cmsearch($opts, $constants, $stats, $gc, $log, "tRNAscan-SE",
                                                               $r_cm_hit_temp, $tmp_trnaseq_file, \$trna_ct, \@rescan_trna_hits);
                            }
                            
                            if (scalar(@rescan_trna_hits) > 0) {
                                if ($rescan_trna_hits[0]->{score} > $best_score) {
                                    $best_score = $rescan_trna_hits[0]->{score};
                                    $total_intron_len += $r_intron_hit->{intron_len};
                                    $r_cm_hit->{acodon} = $rescan_trna_hits[0]->{acodon};
                                    $r_cm_hit->{acodon_pos} = $rescan_trna_hits[0]->{acodon_pos};
                                    $r_cm_hit->{isotype} = $rescan_trna_hits[0]->{isotype};
                                    $r_cm_hit->{score} = $rescan_trna_hits[0]->{score};
                                    $r_cm_hit->{hmm_score} = $rescan_trna_hits[0]->{hmm_score};
                                    $r_cm_hit->{ss_score} = $rescan_trna_hits[0]->{ss_score};
                                    $r_cm_hit->{is_pseudo} = $rescan_trna_hits[0]->{is_pseudo};
                                    $r_cm_hit->{ss} = $rescan_trna_hits[0]->{ss};
                                    $r_cm_hit->{seq} = $rescan_trna_hits[0]->{seq};
#                                    $r_cm_hit->{precursor} = substr($r_cm_hit_temp->{precursor}, $rescan_trna_hits[0]->{start} - 1, $rescan_trna_hits[0]->{len} + $total_intron_len);
                                    $r_cm_hit->{model} = $rescan_trna_hits[0]->{model};
                                    $r_cm_hit->{len} = $rescan_trna_hits[0]->{len} + $total_intron_len;
                                    $new_start = $rescan_trna_hits[0]->{start};
                                    if ($r_intron_hit->{overlap} > -1) {
                                        $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{seq} = $r_intron_hit->{intron_seq};
                                        $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{start} = $r_intron_hit->{start};
                                        $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{end} = $r_intron_hit->{end};
                                        $r_cm_hit->{introns}->[$r_intron_hit->{overlap}]->{type} = "NCI";
                                    }
                                    else {
                                        push(@{$r_cm_hit->{introns}}, {seq=>$r_intron_hit->{intron_seq}, start=>$r_intron_hit->{start},
                                                                       end=>$r_intron_hit->{end}, type=>"NCI"});
                                    }
                                    if (scalar(@{$rescan_trna_hits[0]->{introns}}) > 0) {
                                        $r_new_intron = $rescan_trna_hits[0]->{introns}->[0];
                                        my $set_ci = 0;
                                        foreach $r_intron (@{$r_cm_hit->{introns}}) {
                                            if ($r_intron->{type} eq "CI") {
                                                $r_intron = $r_new_intron;
                                                $set_ci = 1;
                                            }
                                        }
                                        if (!$set_ci) {
                                            push(@{$r_cm_hit->{introns}}, {seq=>$r_new_intron->{seq}, start=>$r_new_intron->{start},
                                                                           end=>$r_new_intron->{end}, type=>"CI"});                                        
                                        }
                                    }
                                    else {
                                        foreach $r_intron (@{$r_cm_hit->{introns}}) {
                                            if ($r_intron->{type} eq "CI") {
                                                $r_intron = undef;
                                            }
                                        }                                    
                                    }
                                }
                            }
                        }
                    }
                }
                if ($total_intron_len > 0) {
                    if ($r_cm_hit->{strand}) {
                        $r_cm_hit->{start} = $r_cm_hit->{start} - length($r_cm_hit->{upstream}) + $new_start - 1;
                        $r_cm_hit->{end} = $r_cm_hit->{start} + $r_cm_hit->{len} - 1;
                    }
                    else {
                        $r_cm_hit->{start} = $r_cm_hit->{start} + length($r_cm_hit->{upstream}) - $new_start + 1;
                        $r_cm_hit->{end} = $r_cm_hit->{start} - $r_cm_hit->{len} + 1;                                        
                    }
                    foreach $r_intron (@{$r_cm_hit->{introns}}) {
                        if (defined $r_intron) {
                            $r_intron->{start} = index($r_cm_hit->{precursor}, uc($r_intron->{seq})) + 1;
                            $r_intron->{end} = $r_intron->{start} + length($r_intron->{seq}) - 1;
                            if (($r_intron->{start} == ($r_cm_hit->{acodon_pos}+4)) && ($r_intron->{type} eq "NCI")) {
                                $r_intron->{type} = "CI";
                            }
                        }
                    }
                    @{$r_cm_hit->{introns}} = sort sort_intron_by_start @{$r_cm_hit->{introns}};
                    $last_end = -1;
                    for (my $ct = 0; $ct < scalar(@{$r_cm_hit->{introns}}); $ct++) {
                        if ($last_end == $r_cm_hit->{introns}->[$ct]->{start} - 1) {
                            $last_end = $r_cm_hit->{introns}->[$ct]->{end};
                            $r_cm_hit->{introns}->[$ct]->{start} = $r_cm_hit->{introns}->[$ct - 1]->{start};
                            $r_cm_hit->{introns}->[$ct]->{seq} = $r_cm_hit->{introns}->[$ct - 1]->{seq} . $r_cm_hit->{introns}->[$ct]->{seq};
                            $r_cm_hit->{introns}->[$ct - 1] = undef;
                        }
                        $last_end = $r_cm_hit->{introns}->[$ct]->{end};
                    }
                }                
            }
        }
    }
    @$r_sec_pass_hits = sort sort_cm_hits_for_output @$r_sec_pass_hits;
}

# Run Infernal cmsearch for noncanonical introns
sub run_cmsearch_intron {
    
    my $self = shift;
    my $opts = shift;
    my $log = shift;
    my ($r_intron_hit_list, $r_cm_file, $r_cm_hit, $tmp_trnaseq_file, $r_over_cutoff_count) = @_;
    
    my ($cms_output, $r_intron_hit, $ct, $intronDesc);
    my ($pre_intron, $intron, $post_intron, $intron_start, $intron_span);

    if (!$self->exec_cmsearch(0, $r_cm_file, 0, $tmp_trnaseq_file, $r_cm_hit->{seqname}, \$cms_output)) {
        return 0;
    }
    $self->parse_cmsearch($cms_output, $r_intron_hit_list, 0, $r_cm_hit, 0);
    $$r_over_cutoff_count = 0;
    foreach $r_intron_hit (@$r_intron_hit_list) {
        $ct++;
        if ($r_intron_hit->{score} >= $self->{BHB_cm_cutoff}) {
            $$r_over_cutoff_count++;
            
            if ($r_intron_hit->{ss} =~ /^([\<\-\.]{11,})(\-\<[<.]+[_.]{4,}[>.]{9,}\-[.]*\-)([-.>]+)$/) {
                $pre_intron  = $1;
                $intron      = $2;
                $post_intron = $3;
                $intron_start = length($pre_intron);
                $intron_span = length($intron);
            }
            $r_intron_hit->{intron_seq} = substr($r_intron_hit->{seq}, $intron_start, $intron_span);                 
            $r_intron_hit->{intron_seq} = uc($r_intron_hit->{intron_seq});
            $r_intron_hit->{intron_seq} =~ s/U/T/g;
            $r_intron_hit->{intron_seq} =~ s/-//g;
            $r_intron_hit->{intron_len} = length($r_intron_hit->{intron_seq});

            $r_intron_hit->{subseq_start} += $intron_start;
            $r_intron_hit->{subseq_end} -= length($post_intron);
            $r_intron_hit->{start} = $r_intron_hit->{subseq_start};
            $r_intron_hit->{end} = $r_intron_hit->{subseq_end};                
            
            if ($r_intron_hit->{strand}) {
                $r_intron_hit->{tRNA_start} += $intron_start;	
                $r_intron_hit->{tRNA_end} -= length($post_intron);
            }
            else {
                $r_intron_hit->{tRNA_start} -= $intron_start;	
                $r_intron_hit->{tRNA_end} += length($post_intron);
            }
            $r_intron_hit->{start} -= length($r_cm_hit->{upstream});
            $r_intron_hit->{end} -= length($r_cm_hit->{upstream});                
            if ($r_intron_hit->{strand}) {
                $r_intron_hit->{tRNA_start} -= length($r_cm_hit->{upstream});	
                $r_intron_hit->{tRNA_end} -= length($r_cm_hit->{upstream});
            }
            else {
                $r_intron_hit->{tRNA_start} += length($r_cm_hit->{upstream});	
                $r_intron_hit->{tRNA_end} += length($r_cm_hit->{upstream});
            }
        }
        else {
            $log->write_line("Low infernal score for noncanonical intron detection of $r_cm_hit->{ID} intron$ct: $r_intron_hit->{score}");
            $log->write_line("CMSearch Hit#$ct: $r_intron_hit->{tRNA_start}-$r_intron_hit->{tRNA_end},".
                " Sc: $r_intron_hit->{score},  Len: ".(abs($r_intron_hit->{tRNA_start} - $r_intron_hit->{tRNA_end}) + 1));
        }
    }        

    return 1;
}

# Run Infernal cmsearch for genome-wide missing tRNA scanning, return results in $r_cms_hit_list array reference

sub run_gw_cmsearch {

    my $self = shift;
    my $opts = shift;
    my $log = shift;
    my ($r_cms_hit_list, $r_cur_cm_file, $tmp_seq_file, $seqname, $halves) = @_;
    
    my ($scan_len, $cms_output, $r_cms_hit, $over_cutoff, $ct, $trnaDesc);
    
    my $cm_file_input = $$r_cur_cm_file;
    $self->set_search_params($opts, \$scan_len, $r_cur_cm_file, $self->{max_cmsearch_tRNA_length},
                     $self->{max_tRNA_length}, "", 0);
    
    # run cmsearch
    if ($halves) {
        $$r_cur_cm_file = $cm_file_input;
    }
    else {
        $$r_cur_cm_file = $self->{arch_gw_scan_cm_file_path};
    }
    if (!$self->exec_cmsearch($scan_len, $r_cur_cm_file, 1, $tmp_seq_file, $seqname, \$cms_output)) {
        return 0;
    }
    $self->parse_gw_cmsearch($cms_output, $r_cms_hit_list, 0);
        
    # Go thru hit list, save info for tRNA hits with sub-cutoff scores

    if (!$halves) {
        $ct = 0;
        $over_cutoff = 0;
    
        foreach $r_cms_hit (@$r_cms_hit_list) {
            $ct++;
            if ($r_cms_hit->{score} >= $self->{cm_cutoff}) {
                $over_cutoff++;
            }
            else {
                $log->write_line("Low cmsearch score for $ct: $r_cms_hit->{score}");
                $trnaDesc = "(CMSearch Hit#$ct: $r_cms_hit->{tRNA_start}-$r_cms_hit->{tRNA_end},".
                    " Sc: $r_cms_hit->{score},  Len: ".(abs($r_cms_hit->{tRNA_start} - $r_cms_hit->{tRNA_end}) + 1).") ";
                if ($opts->save_falsepos()) {
                   &write_tRNA($opts->falsepos_file(), $ct, $trnaDesc, $r_cms_hit->{seq}, 0);
                }           
           }
        }        
    
        # report if no scores over 0 bit reporting threshold
    
        if ($over_cutoff == 0) {
            if (!$opts->results_to_stdout()) {
                $log->write_line("No extra CMSearch hits above cutoff found for $seqname");
            }
        }
        else {
            $log->write_line("Found ".$over_cutoff." extra Infernal hits.");        
        }
    }
    return 1;
}

# Parse the text results of a genome-wide cmsearch into data fields

sub parse_gw_cmsearch {

    my $self = shift;
    my ($cms_output, $r_cms_hit_list, $overlaps_OK) = @_;
    
    my (@cms_lines, $subseq_start, $subseq_end, $hit_overlap,
	   $cur_seq_name, $line_ct, $cms_hit_ct, $score, $strand, $skip,
	   $top_score, $overlapping_hit, $prev_start, $prev_end, $prev_seq_name);
    
    @cms_lines = split(/\n/, $cms_output);
    $cms_hit_ct = -1;
    $score   = 0;

    for ($line_ct = 0; $line_ct <= $#cms_lines; $line_ct++) 
    {
        if ($cms_lines[$line_ct] =~ /^>(\S.+)$/) {  
            $prev_seq_name = $cur_seq_name;
            $cur_seq_name = $1;
            $skip = 0;
            
            if ($cur_seq_name ne $prev_seq_name) {
                $top_score  = -1;
                $prev_start = -1;
                $prev_end   = -1;
            }
        }
        elsif ($cms_lines[$line_ct] =~ /^\s*Plus strand results/)
        {
            $strand = 1;
        }
        elsif ($cms_lines[$line_ct] =~ /^\s*Minus strand results/)
        {
            $strand = 0;
        }
        elsif ($cms_lines[$line_ct] =~ /^\s*Query =\s*(\d+)\s*\-\s*(\d+), Target =\s*(\d+)\s*\-\s*(\d+)/) 
        {
            $skip = 0;
            $subseq_start  = $3;
            $subseq_end    = $4;
            
            $hit_overlap = eval(!$overlaps_OK &&
                    (($cur_seq_name eq $prev_seq_name) && 
                     (seg_overlap($subseq_start, $subseq_end, $prev_start, $prev_end))));
        }
#        elsif ($cms_lines[$line_ct] =~ /^\s*Score =\s*([0-9.\-]+), E =\s*([0-9e.\-]+), P =\s*([0-9e.\-]+), GC =\s*(\d+)/) 
        elsif ($cms_lines[$line_ct] =~ /^\s*Score =\s*([0-9.\-]+), /) 
        {
            $score = $1;
                
            # If true, Don't save current hit -- Advance to next set of alignment info
            if (($score < $top_score) && ($hit_overlap))
            {
                $line_ct += 4;
                $skip = 1;
                next;
            }		
            
            # Save current hit
            # If it's a new sequence, or the same seq but no overlap, save the last hit
            # by advancing the hit counter
            if (!$hit_overlap) 
            {
                # convert cmsearch secondary structure
                if ($cms_hit_ct > -1) {
                    ($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq}) =
                        $self->format_cmsearch_output($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq});
                }

                $cms_hit_ct++;
                push(@$r_cms_hit_list,
                     {hit_seqname => "", score=>-1, ss=>"", seq=>"", model=>"",
                      tRNA_start=>-1, tRNA_end=>-1, 
                      strand=>$strand, hit_source=>"Infernal"});
            }
            
            # Save current hit as best non-overlapping hit
            $prev_start  = $subseq_start;
            $prev_end    = $subseq_end;
            $top_score   = $score;
            
            $r_cms_hit_list->[$cms_hit_ct]->{hit_seqname}= $cur_seq_name;
            $r_cms_hit_list->[$cms_hit_ct]->{score}  = $score;
            $r_cms_hit_list->[$cms_hit_ct]->{ss}     = "";
            $r_cms_hit_list->[$cms_hit_ct]->{seq}    = "";
            $r_cms_hit_list->[$cms_hit_ct]->{model}  = "";
            $r_cms_hit_list->[$cms_hit_ct]->{subseq_start}  = $subseq_start;
            $r_cms_hit_list->[$cms_hit_ct]->{subseq_end}  = $subseq_end;
            $r_cms_hit_list->[$cms_hit_ct]->{tRNA_start}  = $subseq_start;
            $r_cms_hit_list->[$cms_hit_ct]->{tRNA_end}  = $subseq_end;
            $r_cms_hit_list->[$cms_hit_ct]->{tRNA_len} = abs($subseq_end - $subseq_start) + 1;
        }
	
        # Parse model structure line 
	
        elsif ($cms_lines[$line_ct] =~ /^\s+([(),<>._\-,\[\]\{\}\:]{1,60})/)
        {
            if (!$skip) {
                $r_cms_hit_list->[$cms_hit_ct]->{ss} .= $1;
            
                # Parse model sequence line
                if ($cms_lines[$line_ct + 1] =~ /^\s+\d+\s+([a-zA-Z\-]{1,60})\s+\d+/) {  
                    $r_cms_hit_list->[$cms_hit_ct]->{model} .= $1;  
                } 
                # Parse target sequence line
                if ($cms_lines[$line_ct + 3] =~ /^\s+\d+\s+([a-zA-Z\-]{1,60})\s+\d+/) {  
                    $r_cms_hit_list->[$cms_hit_ct]->{seq} .= $1;  
                }
            }
            # Advance to next set of alignment info
            $line_ct += 3;
        }
    }
    # convert cmsearch secondary structure
    if ($cms_hit_ct > -1) {
        ($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq}) =
            $self->format_cmsearch_output($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq});
    }
}

sub parse_covels_hit {

    my $self = shift;
    my ($covels_hit, $r_covels_hit_elements, $ts_start, $sense_strand) = @_;

    my $covels_hit_found = 0;

    if ($covels_hit =~ /^\s*(\S+)\s+(\d+)\s+(\d+).+: (\S+)\s*/o)  {
        $r_covels_hit_elements->{score} = $1;
        $r_covels_hit_elements->{subseq_start} = $2;
        $r_covels_hit_elements->{subseq_end} = $3;
        $r_covels_hit_elements->{hit_seqname} = $4;
        $covels_hit_found = 1;        
    }

    if ($covels_hit_found) {
        if ($sense_strand) {
            $r_covels_hit_elements->{tRNA_start} = $ts_start + $r_covels_hit_elements->{subseq_start} - 1;        
            $r_covels_hit_elements->{tRNA_end} = $ts_start + $r_covels_hit_elements->{subseq_end} -1;
        }
        else {
            $r_covels_hit_elements->{tRNA_start} = $ts_start - $r_covels_hit_elements->{subseq_start} + 1;        
            $r_covels_hit_elements->{tRNA_end} = $ts_start - $r_covels_hit_elements->{subseq_end} + 1;
        }                
        
        return 1;                
    }
    else  {
        return 0;
    }
}                                

# Run covels, return hits in $covels_hit_list array

sub run_covels {

    my $self = shift;
    my $opts = shift;
    my $stats = shift;
    my $log = shift;
    my ($r_covels_hit_list, $r_cur_cm_file, $tmp_trnaseq_file, $r_prescan_tRNA) = @_;

    my ($scan_len, $covels_cmd, $covels_output, $junk, $allhits, $ct,
           $total_hits, $trnaDesc, $report_cutoff, $over_cutoff, $fulltrnaDesc);
    
    my ($covels_hit, $cove_confirmed_ct);
    
    my %covels_hit_elements = ();

    $self->set_search_params($opts, \$scan_len, $r_cur_cm_file, $self->{max_cove_tRNA_length},
                     $r_prescan_tRNA->{len}, $r_prescan_tRNA->{isotype}, 0);
    
    # set covels reporting threshold below 0 (default) if -X param is
    # set below 0 by user

    $report_cutoff = &min(0, $self->{cm_cutoff});
    
    # run Covels

    $covels_cmd = $self->covels_bin()." -w$scan_len -t$report_cutoff $$r_cur_cm_file $tmp_trnaseq_file";
    $covels_output = `$covels_cmd`;

    if (&error_exit_status("Covels-SE", $r_prescan_tRNA->{src_seqname})) {
        print "Exit first loop at 1\n";
        return 0;
    }
    
    ($junk, $allhits) = split(/----------\n\n/, $covels_output);
    @$r_covels_hit_list = split(/\n/, $allhits);

    # count no. of hits over cutoff

    $total_hits = 0;
   
    foreach $covels_hit (@$r_covels_hit_list) {
        %covels_hit_elements = ();
        if (($self->parse_covels_hit($covels_hit, \%covels_hit_elements,
                                     $r_prescan_tRNA->{start}, $r_prescan_tRNA->{strand})) &&
            ($covels_hit_elements{score} >= $self->{cm_cutoff})) {
            $total_hits++;
        }        
    }
    
    # if no tRNAs detected when using a selenocysteine cove model,
    #  try main model and run again before giving up

    if (($total_hits == 0) && 
        (($$r_cur_cm_file eq $self->{Pselc_cm_file_path}) || ($$r_cur_cm_file eq $self->{Eselc_cm_file_path}))) {
        $$r_cur_cm_file = $self->{main_cm_file_path};
        
        # re-run Covels with main model
    
        $covels_cmd = $self->covels_bin()." -w$scan_len -t$report_cutoff $$r_cur_cm_file $tmp_trnaseq_file";
        $covels_output = `$covels_cmd`;
        if (&error_exit_status("Covels-SE", $r_prescan_tRNA->{src_seqname})) {
            print "Exit first loop at 2\n";
            return 0;
        }
        
        ($junk,$allhits) = split(/----------\n\n/,$covels_output);
        @$r_covels_hit_list = split(/\n/, $allhits);
    }

    # Go thru hit list, save info for tRNA hits with sub-cutoff scores

    $ct = 0;
    $over_cutoff = 0;
    $trnaDesc = "";

    foreach $covels_hit (@$r_covels_hit_list) {
        %covels_hit_elements = ();
        if ($self->parse_covels_hit($covels_hit, \%covels_hit_elements,
                                    $r_prescan_tRNA->{start}, $r_prescan_tRNA->{strand})) {
            $ct++;
            if ($covels_hit_elements{score} >= $self->{cm_cutoff}) {
                $over_cutoff++;
            }
            else {
                $log->write_line("Low covels score for $r_prescan_tRNA->{name}.$ct: $covels_hit_elements{score}");
                $trnaDesc .= "(Cove Hit#$ct: $covels_hit_elements{tRNA_start}-$covels_hit_elements{tRNA_end},".
                    " Sc: $covels_hit_elements{score},  Len: ".(abs($covels_hit_elements{tRNA_start} - $covels_hit_elements{tRNA_end}) + 1).") ";
            }
        }
    }        
    
    # report if no scores over 0 bit reporting threshold

    if ($over_cutoff == 0) {
        if ((!$opts->results_to_stdout()) &&
            ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run())) {
            $log->write_line("Covels score(s) below cutoff for $r_prescan_tRNA->{name}. Skipping...");
        }
        if ($opts->save_falsepos()) {
            $fulltrnaDesc = "(Fp Hit: $r_prescan_tRNA->{start}-$r_prescan_tRNA->{end}, ".
                (abs($r_prescan_tRNA->{start} - $r_prescan_tRNA->{end}) + 1)." bp, Src: $r_prescan_tRNA->{hit_source}) ".$trnaDesc;
        
            $stats->increment_fpos_base_ct(length($r_prescan_tRNA->{seq}));          
            &write_tRNA($opts->falsepos_file(), $r_prescan_tRNA->{name}, $fulltrnaDesc, $r_prescan_tRNA->{seq}, 0);
        }           
    }

    return 1;
}

sub run_coves {

    my $self = shift;
    my ($tmp_trnaseq_file, $seq_name, $cm_file) = @_;

    my ($covseq, $covss, $coves_output, $junk, @coves_lines, $sec_struct, $coves_score);
    
    my $coves_cmd = $self->{coves_bin}." -s $cm_file $tmp_trnaseq_file";
    $coves_output = `$coves_cmd`;
    if (&error_exit_status("Coves-SE", $seq_name)) {
        print STDERR "Skipping tRNA anticodon & type prediction\n\n";
        return ("Error", "", -1);
    }

    ($junk, $sec_struct) = split(/----------\n\n/, $coves_output);
    @coves_lines = split(/\n/,$sec_struct);
    $covseq = '';
    $covss = '';
    $coves_score = -1000;
    $seq_name =~ s/(\W)/\\$1/g;

    foreach (@coves_lines) {
        if (/^\s+$seq_name\s([a-zA-Z\-]{1,60})\s*/) {
            $covseq .= $1;  
        } 
        if (/^\s+$seq_name\s([\.\<\>\ ]{1,60})/) {
            $covss .= $1; 
        }
        if (/^\s*(\S+)\sbits\s:\s$seq_name/) {
            $coves_score = $1; 
        }
    }
           
    $covss =~ s/\s//g;     #  take spaces out of alignment        
    $covseq =~ s/-//g;     #  take '-' gaps out of seq

    if (($covseq eq '') || ($covss eq '')) {
        print STDERR "Could not complete coves successfully for $seq_name\n",
        "because unable to parse coves secondary structure string.\n",
        "Skipping tRNA anticodon & type prediction\n";
        return ("Error", "", -1);
    }

    return ($covseq, $covss, $coves_score);
}

sub analyze_with_cove {

    my $self = shift;
    my $opts = shift;
    my $constants = shift;
    my $stats = shift;
    my $gc = shift;
    my $log = shift;
    my $program_id = shift;
    my ($r_prescan_tRNA, $tmp_trnaseq_file, $r_curseq_trnact, $r_sec_pass_hits) = @_;

    my (@covels_hit_list, $cur_cm_file, $covels_hit, $cove_confirmed_ct);
    my ($covseq, $covss, $coves_score, $r_introns);
    my ($cv_anticodon, $acodon_index, $cv_type, $introns, $hmm_score, $ss_score, $pseudo_gene_flag);
    my %covels_hit_elements = ();
    my $cov_hit = {};

    $cove_confirmed_ct = 0;

    if (!$self->run_covels($opts, $stats, $log, \@covels_hit_list, \$cur_cm_file, $tmp_trnaseq_file, $r_prescan_tRNA))  {
        return 0;
    }
        
    # Loop to parse covels tRNA hit(s) and run Coves on each tRNA
    
    foreach $covels_hit (@covels_hit_list) {
        %covels_hit_elements = ();
        if ((!$self->parse_covels_hit($covels_hit, \%covels_hit_elements,
                                      $r_prescan_tRNA->{start}, $r_prescan_tRNA->{strand})) ||
            ($covels_hit_elements{score} < $self->{cm_cutoff})) {
            next; 
        }                       

        $$r_curseq_trnact++;
        
        $covels_hit_elements{upstream} = $r_prescan_tRNA->{upstream};
        $covels_hit_elements{downstream} = $r_prescan_tRNA->{downstream};

        if (($covels_hit_elements{subseq_start} == 1) && ($covels_hit_elements{subseq_end} == $r_prescan_tRNA->{len})) {
            $covels_hit_elements{tRNA_len} = $r_prescan_tRNA->{len};
        }
        else {
            # get correct subseq for coves & save to file
            $covels_hit_elements{tRNA_len} = $covels_hit_elements{subseq_end} - $covels_hit_elements{subseq_start} + 1;
            &write_tRNA($tmp_trnaseq_file, $covels_hit_elements{hit_seqname}, " ",
                        substr($r_prescan_tRNA->{seq}, $covels_hit_elements{subseq_start} - 1, $covels_hit_elements{tRNA_len}), 1);
            if ($covels_hit_elements{subseq_start} > 1) {
                $covels_hit_elements{upstream} .= substr($r_prescan_tRNA->{seq}, 0, $covels_hit_elements{subseq_start} - 1);
            }
            if ($covels_hit_elements{subseq_end} < $r_prescan_tRNA->{len}) {
                $covels_hit_elements{downstream} = substr($r_prescan_tRNA->{seq}, $covels_hit_elements{subseq_end}) .
                    $covels_hit_elements{downstream};
            }
        }                       
        $stats->increment_coves_base_ct($covels_hit_elements{tRNA_len});
    
        $covels_hit_elements{name} = $covels_hit_elements{hit_seqname}.".t".$$r_curseq_trnact;
        
        ($covseq, $covss, $coves_score) = 
            $self->run_coves($tmp_trnaseq_file, $r_prescan_tRNA->{src_seqname}, $cur_cm_file);
        
        # look for intron

        ($cv_anticodon, $acodon_index, $cv_type, $r_introns, $hmm_score, $ss_score, $pseudo_gene_flag) = 
             $self->decode_tRNA_properties ($opts, $gc, $log, $coves_score, $covseq, $covss, $r_prescan_tRNA,
                          $covels_hit_elements{tRNA_start}, $covels_hit_elements{tRNA_end}, $cur_cm_file, $tmp_trnaseq_file);
        
        $cov_hit = {};
        $cov_hit =
             {seqname =>$covels_hit_elements{hit_seqname}, score=>$coves_score, ss=>$covss, seq=>$covseq, model=>"",
              start=>$covels_hit_elements{tRNA_start}, end=>$covels_hit_elements{tRNA_end}, len=>$covels_hit_elements{tRNA_len},
              ID=>$covels_hit_elements{name},
              acodon=>$cv_anticodon, acodon_pos =>$acodon_index, isotype=>$cv_type,
              introns=>$r_introns, hmm_score=>$hmm_score, 
              ss_score=>$ss_score, is_pseudo=>$pseudo_gene_flag,
              src_seqlen=>$r_prescan_tRNA->{src_seqlen}, src_seqname=>$covels_hit_elements{hit_seqname},
              strand=>$r_prescan_tRNA->{strand}, hit_source=>$r_prescan_tRNA->{hit_source},
              upstream=>$covels_hit_elements{upstream}, downstream=>$covels_hit_elements{downstream}, extra=>0};
        
        if (!$self->{CM_check_for_introns}) {
            &output_tRNA($opts, $gc, $log, $self->{tab_results}, $self->{get_hmm_score}, $program_id,
                     $r_prescan_tRNA, $cov_hit, $$r_curseq_trnact);
        
            $cove_confirmed_ct++;
        }
        else {
            push (@$r_sec_pass_hits, $cov_hit);
        }
    }            # while more covels_hits
   
    return $cove_confirmed_ct;
}

# Format command and run Infernal cmsearch

sub exec_cmsearch {

    my $self = shift;
    my ($scan_len, $r_cm_file, $scan_genome, $tmp_trnaseq_file, $seq_name, $r_cms_output) = @_;
    
    my $cm_options = "-g --fil-no-hmm --toponly";
    if ($scan_genome) {
        $cm_options = "-g";
    }
    
    my $cms_cmd = "$self->{cmsearch_bin} $cm_options $$r_cm_file $tmp_trnaseq_file";
    $$r_cms_output = `$cms_cmd`;

    if (&error_exit_status("cmsearch", $seq_name)) {
        print STDERR "Exited at failed cmsearch\n";
        return 0;
    }
    
    return 1;
}

# Run Infernal cmsearch, return results in $r_cms_hit_list array reference

sub run_cmsearch {

    my $self = shift;
    my $opts = shift;
    my $stats = shift;
    my $log = shift;
    my ($r_cms_hit_list, $r_cur_cm_file, $tmp_trnaseq_file, $r_prescan_tRNA) = @_;
    
    my ($scan_len, $cms_cmd, $cms_output, $r_cms_hit, $total_hits, $over_cutoff, $ct, $trnaDesc, $fulltrnaDesc);
    
    $self->set_search_params($opts, \$scan_len, $r_cur_cm_file, $self->{max_cmsearch_tRNA_length},
                     $r_prescan_tRNA->{len}, $r_prescan_tRNA->{isotype}, 0);
    
    # run cmsearch
    
    if (!$self->exec_cmsearch($scan_len, $r_cur_cm_file, 0, $tmp_trnaseq_file, $r_prescan_tRNA->{src_seqname}, \$cms_output)) {
        return 0;
    }
    
    $self->parse_cmsearch($cms_output, $r_cms_hit_list, 0, $r_prescan_tRNA, 1);

    # count no. of hits over cutoff
    
    $total_hits = 0;
    foreach $r_cms_hit (@$r_cms_hit_list) {
        if ($r_cms_hit->{score} >= $self->{cm_cutoff}) {
            $total_hits++;
        }        
    }
    
    # if no tRNAs detected when using a selenocysteine cove model,
    #  try main model and run again before giving up

    if (($total_hits == 0) && 
        (($$r_cur_cm_file eq $self->{Pselc_cm_file_path}) || ($$r_cur_cm_file eq $self->{Eselc_cm_file_path}))) {
        $$r_cur_cm_file = $self->{main_cm_file_path};
        
        # re-run cmsearch with main model
    
        if (!$self->exec_cmsearch($scan_len, $r_cur_cm_file, 0, $tmp_trnaseq_file, $r_prescan_tRNA->{src_seqname}, \$cms_output)) {
            return 0;
        }
        $self->parse_cmsearch($cms_output, $r_cms_hit_list, 0, $r_prescan_tRNA);
    }
    
    # Go thru hit list, save info for tRNA hits with sub-cutoff scores

    $ct = 0;
    $over_cutoff = 0;
    $trnaDesc = "";

    foreach $r_cms_hit (@$r_cms_hit_list) {
        $ct++;
        if ($r_cms_hit->{score} >= $self->{cm_cutoff}) {
            $over_cutoff++;
        }
        else {
            $log->write_line("Low covels score for $r_prescan_tRNA->{name}.$ct: $r_cms_hit->{score}");
            $trnaDesc .= "(CMSearch Hit#$ct: $r_cms_hit->{tRNA_start}-$r_cms_hit->{tRNA_end},".
                " Sc: $r_cms_hit->{score},  Len: ".(abs($r_cms_hit->{tRNA_start} - $r_cms_hit->{tRNA_end}) + 1).") ";
        }
    }        

    # report if no scores over 0 bit reporting threshold

    if ($over_cutoff == 0) {
        if ((!$opts->results_to_stdout()) &&
            ($opts->eufind_mode() || $opts->tscan_mode() || $opts->use_prev_ts_run())) {
            $log->write_line("CMSearch score(s) below cutoff for $r_prescan_tRNA->{name}. Skipping...");
        }
        if ($opts->save_falsepos()) {
            $fulltrnaDesc = "(Fp Hit: $r_prescan_tRNA->{start}-$r_prescan_tRNA->{end}, ".
                (abs($r_prescan_tRNA->{start} - $r_prescan_tRNA->{end}) + 1)." bp, Src: $r_prescan_tRNA->{hit_source}) ".$trnaDesc;
        
            $stats->increment_fpos_base_ct(length($r_prescan_tRNA->{seq}));          
            &write_tRNA($opts->falsepos_file(), $r_prescan_tRNA->{name}, $fulltrnaDesc, $r_prescan_tRNA->{seq}, 0);
        }           
    }

    return 1;
}

# Parse the text results of a cmsearch into data fields

sub parse_cmsearch {

    my $self = shift;
    my ($cms_output, $r_cms_hit_list, $overlaps_OK, $r_prescan_tRNA, $format_struct) = @_;
    
    my (@cms_lines, $subseq_start, $subseq_end, $hit_overlap,
	   $cur_seq_name, $line_ct, $cms_hit_ct, $score, $skip,
	   $top_score, $overlapping_hit, $prev_start, $prev_end, $prev_seq_name);
    
    @cms_lines = split(/\n/, $cms_output);
    $cms_hit_ct = -1;
    $score   = 0;

    for ($line_ct = 0; $line_ct <= $#cms_lines; $line_ct++) 
    {
        if ($cms_lines[$line_ct] =~ /^>(\S.+)$/) {  
            $prev_seq_name = $cur_seq_name;
            $cur_seq_name = $1;
            $skip = 0;
            
            if ($cur_seq_name ne $prev_seq_name) {
                $top_score  = -1;
                $prev_start = -1;
                $prev_end   = -1;
            }
        }
        elsif ($cms_lines[$line_ct] =~ /^\s*Query =\s*(\d+)\s*\-\s*(\d+), Target =\s*(\d+)\s*\-\s*(\d+)/) 
        {
            $skip = 0;
            $subseq_start  = $3;
            $subseq_end    = $4;
            
            $hit_overlap = eval(!$overlaps_OK &&
                    (($cur_seq_name eq $prev_seq_name) && 
                     (seg_overlap($subseq_start, $subseq_end, $prev_start, $prev_end))));
        }
#        elsif ($cms_lines[$line_ct] =~ /^\s*Score =\s*([0-9.\-]+), E =\s*([0-9e.\-]+), P =\s*([0-9e.\-]+), GC =\s*(\d+)/) 
        elsif ($cms_lines[$line_ct] =~ /^\s*Score =\s*([0-9.\-]+), /) 
        {
            $score = $1;
                
            # If true, Don't save current hit -- Advance to next set of alignment info
            if (($score < $top_score) && ($hit_overlap))
            {
                $line_ct += 4;
                $skip = 1;
                next;
            }		
            
            # Save current hit
            # If it's a new sequence, or the same seq but no overlap, save the last hit
            # by advancing the hit counter
            if (!$hit_overlap) 
            {
                # convert cmsearch secondary structure
                if (($cms_hit_ct > -1) && $format_struct) {
                    ($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq}) =
                        $self->format_cmsearch_output($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq});
                }
                
                $cms_hit_ct++;
                push(@$r_cms_hit_list,
                     {hit_seqname => "", score=>-1, ss=>"", seq=>"", model=>"",
                      tRNA_start=>-1, tRNA_end=>-1, 
                      strand=>$r_prescan_tRNA->{strand}, hit_source=>$r_prescan_tRNA->{hit_source}});
            }
            
            # Save current hit as best non-overlapping hit
            $prev_start  = $subseq_start;
            $prev_end    = $subseq_end;
            $top_score   = $score;
            
            $r_cms_hit_list->[$cms_hit_ct]->{hit_seqname}= $cur_seq_name;
            $r_cms_hit_list->[$cms_hit_ct]->{score}  = $score;
            $r_cms_hit_list->[$cms_hit_ct]->{ss}     = "";
            $r_cms_hit_list->[$cms_hit_ct]->{seq}    = "";
            $r_cms_hit_list->[$cms_hit_ct]->{model}  = "";
            $r_cms_hit_list->[$cms_hit_ct]->{subseq_start}  = $subseq_start;
            $r_cms_hit_list->[$cms_hit_ct]->{subseq_end}  = $subseq_end;
            
            if ($r_prescan_tRNA->{strand}) {
                $r_cms_hit_list->[$cms_hit_ct]->{tRNA_start} = $r_prescan_tRNA->{start} + $subseq_start - 1;	
                $r_cms_hit_list->[$cms_hit_ct]->{tRNA_end}   = $r_prescan_tRNA->{start} + $subseq_end - 1;
            }
            else {
                $r_cms_hit_list->[$cms_hit_ct]->{tRNA_start} = $r_prescan_tRNA->{start} - $subseq_start + 1;	
                $r_cms_hit_list->[$cms_hit_ct]->{tRNA_end}   = $r_prescan_tRNA->{start} - $subseq_end + 1;
            }		
        }
	
        # Parse model structure line 
	
        elsif ($cms_lines[$line_ct] =~ /^\s+([(),<>._\-,\[\]\{\}\:]{1,60})/)
        {
            if (!$skip) {
                $r_cms_hit_list->[$cms_hit_ct]->{ss} .= $1;
            
                # Parse model sequence line
                if ($cms_lines[$line_ct + 1] =~ /^\s+\d+\s+([a-zA-Z\-]{1,60})\s+\d+/) {  
                    $r_cms_hit_list->[$cms_hit_ct]->{model} .= $1;  
                } 
                # Parse target sequence line
                if ($cms_lines[$line_ct + 3] =~ /^\s+\d+\s+([a-zA-Z\-]{1,60})\s+\d+/) {  
                    $r_cms_hit_list->[$cms_hit_ct]->{seq} .= $1;  
                }
            }
             # Advance to next set of alignment info
            $line_ct += 3;
        }
    }
    # convert cmsearch secondary structure
    if (($cms_hit_ct > -1) && $format_struct) {
        ($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq}) =
            $self->format_cmsearch_output($r_cms_hit_list->[$cms_hit_ct]->{ss}, $r_cms_hit_list->[$cms_hit_ct]->{seq});
    }

}

sub format_cmsearch_output {
    
    my $self = shift;
    my $cmsearch_ss = shift;
    my $cmsearch_seq = shift;
    
    $cmsearch_seq =~ s/U/T/g; 
    $cmsearch_seq =~ s/u/t/g;
    for (my $index = 0; $index < length($cmsearch_seq); $index++) {
        if (substr($cmsearch_seq, $index, 1) eq '-') {
            substr($cmsearch_seq, $index, 1) = '*';
            if (length($cmsearch_ss) > $index) {
                substr($cmsearch_ss, $index, 1) = '*';
            }
        }
    }
    $cmsearch_seq =~ s/\*//g;
    $cmsearch_ss =~ s/\*//g;

    $cmsearch_ss =~ s/[,_\-:]/./g;
    $cmsearch_ss =~ s/[>)]/@/g;
    $cmsearch_ss =~ s/[(<]/>/g;
    $cmsearch_ss =~ s/@/</g;
    
    my $diff = length($cmsearch_seq) - length($cmsearch_ss);
    for (my $ct = 0; $ct < $diff; $ct++) {
        $cmsearch_ss .= ".";
    }
    
    return ($cmsearch_ss, $cmsearch_seq);
}

# Runs Infernal cmsearch, and returns the highest score

sub cmsearch_bestscore {

    my $self = shift;
    my $opts = shift;
    my ($tmp_trnaseq_file, $seq_name, $tRNA_len, $cm_file) = @_;

    my ($cms_output, @cms_hit_list, $scan_len, $cur_cm_file, @cms_lines,
        $subseq_start, $subseq_end, $score, $besthit_score, $line,
	$besthit_start, $besthit_end, $hit_ct);

    $self->set_search_params($opts, \$scan_len, \$cur_cm_file, $self->{max_cmsearch_tRNA_length}, $tRNA_len, '', 1);

    if (!$self->exec_cmsearch($scan_len, \$cur_cm_file, 0, $tmp_trnaseq_file, $seq_name, \$cms_output)) {
        return 0;
    }

    @cms_lines = split(/\n/, $cms_output);
    $hit_ct = 0;
    $score = 0;
    $besthit_score = 0;

    foreach $line (@cms_lines)
    {
	if ($line =~ /^\s*Query =\s+(\d+)\s+-\s+(\d+), Target =\s+(\d+)\s+-\s+(\d+)/) 
	{
	    $subseq_start  = $3;
	    $subseq_end    = $4;
	    $hit_ct++;
        }
        elsif ($line =~ /^\s*Score =\s+([0-9.\-]+), GC =\s+(\d+)/) 
	{           
	    $score = $1;
	    if ($score > $besthit_score) {
		$besthit_score  = $score;   
		$besthit_start  = $subseq_start;
		$besthit_end    = $subseq_end;
	    }
	}
    }
    return ($besthit_score, $besthit_start, $besthit_end, $hit_ct); 
}

sub analyze_with_cmsearch {
 
    my $self = shift;
    my $opts = shift;
    my $constants = shift;
    my $stats = shift;
    my $gc = shift;
    my $log = shift;
    my $program_id = shift;
    my ($r_prescan_tRNA, $tmp_trnaseq_file, $r_curseq_trnact, $r_sec_pass_hits) = @_;

    my ($cms_confirmed_ct, $cms_output, $cur_cm_file, @cms_hit_list, $cms_hit, $ct);
    my ($cm_anticodon, $acodon_index, $cm_type, $r_introns, $hmm_score, $ss_score, $pseudo_gene_flag);
    
    my $cm_hit = {};
    
    $cms_confirmed_ct = 0;
    
    if (!$self->run_cmsearch($opts, $stats, $log, \@cms_hit_list, \$cur_cm_file, $tmp_trnaseq_file, $r_prescan_tRNA)) {
        return 0;
    }

    # Loop to process each cmsearch tRNA hit
    
    foreach $cms_hit (@cms_hit_list) {
        
        if ($cms_hit->{score} < $self->{cm_cutoff}) {
            next;
        }
        
        $$r_curseq_trnact++;

        $cms_hit->{upstream} = $r_prescan_tRNA->{upstream};
        $cms_hit->{downstream} = $r_prescan_tRNA->{downstream};

        if (($cms_hit->{subseq_start} == 1) && ($cms_hit->{subseq_end} == $r_prescan_tRNA->{len})) {
            $cms_hit->{tRNA_len} = $r_prescan_tRNA->{len};
        }
        else {
            # get correct subseq & save to file
            if ((length($r_prescan_tRNA->{seq}) >= $cms_hit->{subseq_end} + 2) &&
                (uc(substr($r_prescan_tRNA->{seq}, $cms_hit->{subseq_end}, 3)) eq "CCA") &&
                (substr($cms_hit->{ss}, length($cms_hit->{ss}) - 4) ne "....")) {
                $cms_hit->{subseq_end} += 3;
                if ($cms_hit->{strand}) {
                    $cms_hit->{tRNA_end} += 3;
                }
                else {
                    $cms_hit->{tRNA_end} -= 3;
                }
                $cms_hit->{seq} .= "CCA";
            }
            $cms_hit->{tRNA_len} = $cms_hit->{subseq_end} - $cms_hit->{subseq_start} + 1;
            &write_tRNA($tmp_trnaseq_file, $cms_hit->{hit_seqname}, " ",
                        substr($r_prescan_tRNA->{seq}, $cms_hit->{subseq_start} - 1, $cms_hit->{tRNA_len}), 1);

            if (uc(substr($r_prescan_tRNA->{seq}, $cms_hit->{subseq_start} - 1, $cms_hit->{tRNA_len})) ne uc($cms_hit->{seq})) {
                $cms_hit->{seq} = substr($r_prescan_tRNA->{seq}, $cms_hit->{subseq_start} - 1, $cms_hit->{tRNA_len});
            }
            if ($cms_hit->{subseq_start} > 1) {
                $cms_hit->{upstream} .= substr($r_prescan_tRNA->{seq}, 0, $cms_hit->{subseq_start} - 1);
            }
            if ($cms_hit->{subseq_end} < $r_prescan_tRNA->{len}) {
                $cms_hit->{downstream} = substr($r_prescan_tRNA->{seq}, $cms_hit->{subseq_end}) . $cms_hit->{downstream};
            }
        }                       

        $cms_hit->{name}  = $cms_hit->{hit_seqname}.".t$$r_curseq_trnact";
        
        ($cm_anticodon, $acodon_index, $cm_type, $r_introns, $hmm_score, $ss_score, $pseudo_gene_flag) = 
             $self->decode_tRNA_properties ($opts, $gc, $log, $cms_hit->{score}, $cms_hit->{seq}, $cms_hit->{ss}, $r_prescan_tRNA,
                          $cms_hit->{tRNA_start}, $cms_hit->{tRNA_end}, $cur_cm_file, $tmp_trnaseq_file);
        
        $cm_hit = {};
        $cm_hit =
             {seqname =>$cms_hit->{hit_seqname}, score=>$cms_hit->{score}, ss=>$cms_hit->{ss}, seq=>$cms_hit->{seq}, model=>$cms_hit->{model},
              start=>$cms_hit->{tRNA_start}, end=>$cms_hit->{tRNA_end}, len=>$cms_hit->{tRNA_len}, ID=>$cms_hit->{name},
              acodon=>$cm_anticodon, acodon_pos =>$acodon_index, isotype=>$cm_type,
              introns=>$r_introns, hmm_score=>$hmm_score, 
              ss_score=>$ss_score, is_pseudo=>$pseudo_gene_flag,
              src_seqlen=>$r_prescan_tRNA->{src_seqlen}, src_seqname=>$cms_hit->{hit_seqname},
              strand=>$r_prescan_tRNA->{strand}, hit_source=>$r_prescan_tRNA->{hit_source}, 
              upstream=>$cms_hit->{upstream}, downstream=>$cms_hit->{downstream}, extra=>0};
        
        if (!$self->{CM_check_for_introns}) {
            &output_tRNA($opts, $gc, $log, $self->{tab_results}, $self->{get_hmm_score}, $program_id,
                         $r_prescan_tRNA, $cm_hit, $$r_curseq_trnact);	
        
            $cms_confirmed_ct++;
        }
        else {
            push(@$r_sec_pass_hits, $cm_hit);
        }
    }	# while more cmsearch hits
    
    return $cms_confirmed_ct;
}	    

sub sort_cm_hits_by_start {
    
    my $self = shift;
    my $cms_hits = shift;
    
    my $a_start = $a->{start};
    my $b_start = $b->{start};
    
    if ($a->{strand} == 0) {
        $a_start = $a->{end};
    }
    if ($b->{strand} == 0) {
        $b_start = $b->{end};
    }
    
    return ($a->{seqname} cmp $b->{seqname} ||
            $a_start <=> $b_start);    
}

sub sort_cm_hits_for_output {
    
    my $self = shift;
    my $cms_hits = shift;
    
    if ((($a->{strand} == $b->{strand}) && ($a->{strand} == 1)) ||
        ($a->{strand} != $b->{strand})) {
        return ($a->{seqname} cmp $b->{seqname} ||
                $b->{strand} <=> $a->{strand} ||
                $a->{start} <=> $b->{start});
    }
    if (($a->{strand} == $b->{strand}) && ($a->{strand} == 0)) {
        return ($a->{seqname} cmp $b->{seqname} ||
                $b->{strand} <=> $a->{strand} ||
                $b->{start} <=> $a->{start});        
    }
}

sub sort_intron_by_start {
    
    my $self = shift;
    my $introns = shift;
        
    return ($a->{start} <=> $b->{end});    
}

1;
