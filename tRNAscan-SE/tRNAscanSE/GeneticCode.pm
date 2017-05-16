# tRNAscanSE/GeneticCode.pm
# This class describes the genetic codes used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::GeneticCode;

use strict;
use tRNAscanSE::Utils;

sub new {
    my $class = shift;
    my $self = {};

    initialize($self);

    bless ($self, $class);
    return $self;
}

sub DESTROY {
    my $self = shift;
}

sub initialize {
    my $self = shift;

    $self->{undef_anticodon} = "???";
    $self->{undef_isotype}   = "Undet";

    my @isotypes = ('Ala', 'Gly', 'Pro', 'Thr', 'Val', 
                 'Ser', 'Arg', 'Leu',
                 'Phe','Asn', 'Lys', 'Asp', 'Glu', 'His', 'Gln', 
                 'Ile', 'Met', 'Tyr', 'Supres', 'Cys', 'Trp', 'SelCys');
    $self->{isotypes} = \@isotypes;
    
    # Amino acid -> Anti-codon list for printing out global tRNA summary

    my %ac_list = (
               'Ala' => [qw/AGC GGC CGC TGC/],
               'Gly' => [qw/ACC GCC CCC TCC/],
               'Pro' => [qw/AGG GGG CGG TGG/],
               'Thr' => [qw/AGT GGT CGT TGT/],
               'Val' => [qw/AAC GAC CAC TAC/],
               
               'Ser' => [qw/AGA GGA CGA TGA ACT GCT/],
               'Arg' => [qw/ACG GCG CCG TCG CCT TCT/],
               'Leu' => [qw/AAG GAG CAG TAG CAA TAA/],
               
               'Phe' => [qw/AAA GAA &nbsp &nbsp /],
               
               'Asn' => [qw/ATT GTT &nbsp &nbsp /],
               'Lys' => [qw/&nbsp &nbsp CTT TTT/],
               
               'Asp' => [qw/ATC GTC &nbsp &nbsp /],
               'Glu' => [qw/&nbsp &nbsp CTC TTC/],
               
               'His' => [qw/ATG GTG &nbsp &nbsp /],
               'Gln' => [qw/&nbsp &nbsp CTG TTG/],
               
               'Tyr' => [qw/ATA GTA &nbsp &nbsp /],
               'Supres' => [qw/&nbsp &nbsp CTA TTA/],
               
               'Ile' => [qw/AAT GAT &nbsp TAT/],
               'Met' => [qw/&nbsp &nbsp CAT &nbsp/],
               
               'Cys' => [qw/ACA GCA &nbsp &nbsp /],
               'Trp' => [qw/&nbsp &nbsp CCA &nbsp/],
               'SelCys' => [qw/&nbsp &nbsp &nbsp TCA/]
               );
    $self->{ac_list} = \%ac_list;    

    $self->{trans_map} = +{};
    $self->{one_let_trans_map} = +{};
}

sub undef_anticodon
{
    my $self = shift;
    return $self->{undef_anticodon};
}

sub undef_isotype
{
    my $self = shift;
    return $self->{undef_isotype};
}

sub isotypes
{
    my $self = shift;
    return $self->{isotypes};
}

sub ac_list
{
    my $self = shift;
    return $self->{ac_list};
}

sub one_let_trans_map
{
    my $self = shift;
    return $self->{one_let_trans_map};
}

sub read_transl_table {
    
    my $self = shift;
    my $opts = shift;
    my $alt_gcode = $opts->alt_gcode();
    my $gc_file = $opts->gc_file();
    
    my %ambig_trans_map = ();
    my %alt_trans_map = ();
    my ($acodon, @expanded_set, $expanded_ac, $gc_file_path);
    
    # Read in default genetic code table (may contain ambiguous bases) at
    # end of this source file

    while (<DATA>) {                
        if ((/^[^\#]/) && 
            (/^([ACGTUNRYSWMKBDHV]{3,3})\s+(\S+)\s+(\S)/i)) {
            $acodon = uc($1);
            $ambig_trans_map{&rev_comp_seq($acodon)} = $2;
            $self->{one_let_trans_map}->{$2} = $3;
        } 
    }                

    $self->{one_let_trans_map}->{$self->{undef_isotype}} = "?";
    $self->{one_let_trans_map}->{"SeC(p)"} = "Z";
    $self->{one_let_trans_map}->{"SeC(e)"} = "Z";

    # Convert any ambiguous bases to make all non-ambigous codons
    #  and save translated amino acid

    @expanded_set = ();
    foreach $acodon (sort keys(%ambig_trans_map)) {
        push(@expanded_set, &expand_ambig($acodon));
        foreach $expanded_ac (@expanded_set) {
            $self->{trans_map}->{$expanded_ac} =  $ambig_trans_map{$acodon};  
        }            
        @expanded_set = ();
    }

    if ($alt_gcode) {
    
        if (-r $gc_file) {
            $gc_file_path = $gc_file;
        }
        elsif (-r "/usr/local/lib/tRNAscanSE/".$gc_file) {
            $gc_file_path = "/usr/local/lib/tRNAscanSE/".$gc_file; 
        }
        else {
            die "FATAL: Could not find $gc_file translation codon file\n\n";
        }
    
        open (GC_TABLE, "$gc_file_path") || 
            die "FATAL: Could not find $gc_file translation codon file\n\n";

        # Read in genetic code table (may contain ambiguous bases)
    
        while (<GC_TABLE>) {                
            if ((/^[^\#]/) 
                && (/^([ACGTUNRYSWMKBDHV]{3,3})\s+(\S+)\s+(\S)/i)) {
                $acodon = uc($1);
                $alt_trans_map{&rev_comp_seq($acodon)} = $2;  
                $self->{one_let_trans_map}->{$2} = $3;  
            } 
        }
            close GC_TABLE;
                                   
        # Convert any ambiguous bases to make all non-ambigous codons
        #  and save translated amino acid
    
        @expanded_set = ();
        foreach $acodon (sort keys(%alt_trans_map)) {
            push(@expanded_set, &expand_ambig($acodon));
            foreach $expanded_ac (@expanded_set) {
                $self->{trans_map}->{$expanded_ac} =  $alt_trans_map{$acodon};  
            }            
            @expanded_set = ();
        }
    }    
}

sub get_tRNA_type {

    my $self = shift;
    my $cm = shift;
    my $ac = shift;                         # anticodon to be decoded
    my $cm_file = shift;

    my $Pselc_cm_file_path = $cm->Pselc_cm_file_path();
    my $Eselc_cm_file_path = $cm->Eselc_cm_file_path();
    
    my ($prev_type,$type);

    if ($ac eq $self->{undef_anticodon}) {
        return $self->{undef_isotype};
    }
    elsif ($cm_file eq $Pselc_cm_file_path) {
        return 'SeC(p)';
    }
    elsif ($cm_file eq $Eselc_cm_file_path) {
        return 'SeC(e)';
    }
    else {
        $prev_type = 'INIT';
        foreach my $exp_codon (&expand_ambig($ac)) {
            $type = $self->{trans_map}->{$exp_codon};
            if (($type ne $prev_type) && ($prev_type ne 'INIT')) {
                return $self->{undef_isotype};
            }
            $prev_type = $type;
        }
        return $type;
    }
}

sub expand_ambig {
    
    my ($ac) = @_;

    $ac = " ".$ac." ";
    
    while (index($ac, 'N') != -1) {
        $ac =~ s/(.*)\s(\S*)N(\S*)\s(.*)/$1 $2A$3 $2C$3 $2G$3 $2T$3 $4/g;
    }
    &expand2(\$ac, 'Y', 'C', 'T'); &expand2(\$ac, 'R', 'A', 'G'); 
    &expand2(\$ac, 'W', 'A', 'T'); &expand2(\$ac, 'S', 'C', 'G'); 
    &expand2(\$ac, 'M', 'A', 'C'); &expand2(\$ac, 'K', 'G', 'T');
    
    &expand3(\$ac, 'V', 'A', 'C', 'G'); &expand3(\$ac, 'B', 'C', 'G', 'T'); 
    &expand3(\$ac, 'H', 'A', 'C', 'T'); &expand3(\$ac, 'D', 'A', 'G', 'T'); 
    
    $ac = substr($ac, 1);
    return (split(/ /, $ac));
}

sub expand2 {
    
    my ($acodon, $ambig_base, $sub1, $sub2) = @_;
    
    while (index($$acodon, $ambig_base) != -1) {
        $$acodon =~ s/(.*)\s(\S*)$ambig_base(\S*)\s(.*)/$1 $2$sub1$3 $2$sub2$3 $4/g;
    }
}

sub expand3 {
    
    my($acodon, $ambig_base, $sub1, $sub2, $sub3) = @_;

    while (index($$acodon, $ambig_base) != -1) {
        $$acodon =~ s/(.*)\s(\S*)$ambig_base(\S*)\s(.*)/$1 $2$sub1$3 $2$sub2$3 $2$sub3$3 $4/g;
    }
}

1;

__DATA__
GCN        Ala        A
TGY        Cys        C
GAY        Asp        D
GAR        Glu        E
TTY        Phe        F
GGN        Gly        G
CAY        His        H
ATH        Ile        I
AAR        Lys        K
TTR        Leu        L
CTN        Leu        L
ATG        Met        M
AAY        Asn        N
CCN        Pro        P
CAR        Gln        Q
AGR        Arg        R
CGN        Arg        R
AGY        Ser        S
TCN        Ser        S
ACN        Thr        T
GTN        Val        V
TGG        Trp        W
TAY        Tyr        Y
TAR        Sup        ?
TGA        SeC        Z

