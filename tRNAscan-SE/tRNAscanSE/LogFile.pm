# tRNAscanSE/LogFile.pm
# This class defines the log file used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::LogFile;

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
    $self->{file_name} = "";             # name of log file
    $self->{FILE_H} = undef;             # file handle
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub open_file
{
    my $self = shift;
    my $file = shift;
    
    my $success = 0;
    
    if (($file eq "-") || ($file eq "/dev/null"))
    {
        $success = open($self->{FILE_H}, ">$file");
        $self->{file_name} = $file;
    }
    else
    {
        &open_for_write(\$self->{FILE_H}, $file);
        select($self->{FILE_H});
        $|=1;
        $self->{file_name} = $file;
        $success = 1;
    }
    return $success;
}

sub close_file
{
    my $self = shift;
    
    if (defined $self->{FILE_H})
    {
        close($self->{FILE_H});
    }
}

sub write_line
{
    my $self = shift;
    my $line = shift;
    
    my $fh = $self->{FILE_H};
    
    print $fh $line . "\n";
}

1;
