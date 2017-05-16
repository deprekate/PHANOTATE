# tRNAscanSE/Constants.pm
# This class defines global constants used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Constants;

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

    $self->{REALLY_BIG_NUMBER} = 1000000000;   # largest sequence length imaginable

    # Source of first-pass hits table
    # C = Cove, T = tRNAscan, E = EufindtRNA, B = both

    my @source_tab = ('Cv', 'Ts', 'Eu', 'Bo');
    $self->{source_tab} = \@source_tab;
    
    $self->{upstream_len} = 60;
    $self->{downstream_len} = 60;
}

sub REALLY_BIG_NUMBER
{
    my $self = shift;
    return $self->{REALLY_BIG_NUMBER};
}

sub source_tab
{
    my $self = shift;
    return $self->{source_tab};
}

sub upstream_len
{
    my $self = shift;
    return $self->{upstream_len};
}

sub downstream_len
{
    my $self = shift;
    return $self->{downstream_len};
}

sub set_temp_file_names
{
    my $self = shift;
    my $temp_dir = shift;

    $self->{tmp_raw} = &tempname($temp_dir, ".raw");               # for raw tscan output
    $self->{tmp_fa} = &tempname($temp_dir, ".fa");                 # for current fasta seq file
    $self->{tmp_trnaseq_file} = &tempname($temp_dir, ".trna");     # for current tRNA seq 
    $self->{tmp_masked_fa} = &tempname($temp_dir, ".masked.fa");   # for current tRNA seq 
}

sub tmp_raw
{
    my $self = shift;
    return $self->{tmp_raw};
}

sub tmp_fa
{
    my $self = shift;
    return $self->{tmp_fa};
}

sub tmp_trnaseq_file
{
    my $self = shift;
    return $self->{tmp_trnaseq_file};
}

sub tmp_masked_fa
{
    my $self = shift;
    return $self->{tmp_masked_fa};
}

1;
