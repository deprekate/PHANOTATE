#! /usr/bin/perl

# Usage: fasta2gsi.perl <seqfile>
# Creates seqfile.gsi
#
# Create a .gsi sequence database index file.
#
# GSI allows multiple files per index, but fasta2gsi.pl
# creates a GSI index for a single FASTA file.
#
# Part of the SQUID sequence analysis library.
# Copyright (C) 1992-1996 Sean R. Eddy


$seqfile = shift;
$gsifile = $seqfile.".gsi";
$tmpfile = $seqfile.".tmpgsi";

# Library of Perl functions for creating GSI index files.
#
# GSI definition: 
#    1 + <nfiles> + <nkeys> total records.
#    Each record = 38 bytes.
#
#  one header record     :  <"GSI"    (32)> <nfiles (2)> <nkeys (4)> 
#  <nfiles> file records :  <filename (32)> <fileno (2)> <fmt   (4)> 
#  <nkeys>  key records  :  <key      (32)> <fileno (2)> <offset(4)> 
#
# Part of the SQUID sequence analysis library.
# Copyright (C) 1992-1996 Sean R. Eddy


# The following numbers MUST match their counterparts in squid.h
#
$sqd_fmt_genbank = 2;
$sqd_fmt_embl    = 4;
$sqd_fmt_fasta   = 7;
$sqd_fmt_pir     = 12;
 
# Function: GSI_WriteHeader(GSIFILE, $filenum, $keynum)
# 
# Write the header of an open GSI file.
#
sub GSI_WriteHeader {
    local(*GSIFILE, $filenum, $keynum) = @_;
    local($header); 
    $header = pack("a32 n N", "GSI", $filenum, $keynum);
    print GSIFILE $header;
    1;
}

# Function: GSI_WriteFileRecord(GSIFILE, $filename, $idx, $fmt)
#
# Write a file record to an open GSI file.
#
sub GSI_WriteFileRecord {
    local(*GSIFILE, $filename, $idx, $fmt) = @_;
    local($record);
    $record = pack("a32 n N", $filename, $idx, $fmt);
    print GSIFILE $record;
    1;
}
    
# Function: GSI_WriteKeyRecord(GSIFILE, $key, $filenum, $offset)
#
# Write a key record to an open GSI file.
#
sub GSI_WriteKeyRecord {
    local(*GSIFILE, $key, $filenum, $offset) = @_;
    local($record);
    $record = pack("a32 n N", $key, $filenum, $offset);
    print GSIFILE $record;
    1;
}


 

# First pass. Create an unsorted flat text file.
#
$curr_offset = 0;
$recnum      = 0;
print "Calculating offsets for $seqfile...\n";
open(TMPFILE,">$tmpfile");
open(SEQFILE,$seqfile);
while (<SEQFILE>)
{
    if (($key) = /^>\s*(\S+)/)
    {
	print TMPFILE "$key 1 $curr_offset\n";
	$recnum++;
    }
    $curr_offset = tell;
}
close(SEQFILE);
close(TMPFILE);

# Sort the temporary file alphabetically on the key.
print "Sorting the intermediate index file...\n";
system("sort -o $tmpfile $tmpfile");

# Second pass. Convert flat text file to binary GSI.
#
print "Writing the final binary GSI file...\n";
open(GSIFILE,">$gsifile");
&GSI_WriteHeader(GSIFILE, 1, $recnum);
&GSI_WriteFileRecord(GSIFILE, $seqfile, 1, $sqd_fmt_fasta);

open(TMPFILE,$tmpfile);
while (<TMPFILE>)
{
    ($key, $filenum, $offset) = split;
    &GSI_WriteKeyRecord(GSIFILE, $key, $filenum, $offset);
}			
close(TMPFILE);
close(GSIFILE);
unlink $tmpfile;

print "Complete.\n";
print "$gsifile indexes $recnum sequence names.\n";
