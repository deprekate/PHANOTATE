# tRNAscanSE/Sequence.pm
# This class describes a sequence and provides functions for handling fasta files in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#
# Perl code for reading FASTA-formatted sequence files
# SRE, Sat Feb 19 19:10:43 1994

# These subroutines read a FASTA formatted file one sequence at a time.
# Open(filename, open_mode) opens a file for reading or wrting.
# Close() closes it when you're done.
#
# read_fasta() returns 1 on success and 0 on failure (end of file).
# When it returns success, the following variables are set:
#
#       $seq_name        = name of sequence (1st word on FASTA title line)
#       $seq_description = description      (remainder of FASTA title line)
#       $seq_length      = length of sequence
#       $sequence        = sequence, gaps and newlines removed
#
# Modified by TMJL  11/95 for use in tRNAscan-SE

package tRNAscanSE::Sequence;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::Constants;
use tRNAscanSE::Options;


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
    
    $self->{max_seq_buffer} = 1000000;     # Max size of seq buffer read in at once
    $self->{seq_buf_overlap} = 200;        # Nucleotides of overlap between buffers
    $self->{seq_index_inc} = 100000;
    
    $self->{saved_line} = "";
    $self->{buffer_overlap_seq} = "";
    $self->{buffer_end_index} = 0;
    $self->{seq_buf_overrun} = 0;
    $self->{buffer_length} = 0;    
    $self->{key_found} = 0;
    $self->{all_seq_indices} = +[];        # Keeps track of indexing into seqs for fast retreival
        
    $self->{seq_id} = 0;
    $self->{seq_name} = "";
    $self->{seq_description} = "";
    $self->{seq_length} = 0;
    $self->{sequence} = undef;
    
    $self->{seq_name_map} = {};
}

sub file_name
{
    my $self = shift;
    if (@_) { $self->{file_name} = shift; }
    return $self->{file_name};
}

sub key_found
{
    my $self = shift;
    if (@_) { $self->{key_found} = shift; }
    return $self->{key_found};
}

sub seq_id
{
    my $self = shift;
    if (@_) { $self->{seq_id} = shift; }
    return $self->{seq_id};
}

sub seq_name
{
    my $self = shift;
    if (@_) { $self->{seq_name} = shift; }
    return $self->{seq_name};
}

sub get_seq_id_from_name
{
    my $self = shift;
    my $name = shift;
    my $id = -1;
    if (defined $self->{seq_name_map}->{$name})
    {
        $id = $self->{seq_name_map}->{$name};
    }
    return $id;
}

sub seq_description
{
    my $self = shift;
    if (@_) { $self->{seq_description} = shift; }
    return $self->{seq_description};
}

sub seq_length
{
    my $self = shift;
    if (@_) { $self->{seq_length} = shift; }
    return $self->{seq_length};
}

sub sequence
{
    my $self = shift;
    if (@_) { $self->{sequence} = shift; }
    return $self->{sequence};
}

sub release_memory
{
    my $self = shift;
    undef($self->{sequence});
}

sub set_seq_info
{
    my $self = shift;
    $self->{seq_name} = shift;
    $self->{seq_description} = shift;
    $self->{seq_length} = shift;
    $self->{sequence} = shift;
}

sub reset_buffer_ct
{
    my $self = shift;
    $self->{buffer_overlap_seq} = "";
    $self->{buffer_end_index} = 0;
    $self->{seq_buf_overrun} = 0;
}

sub buffer_end_index
{
    my $self = shift;
    if (@_) { $self->{buffer_end_index} = shift; }
    return $self->{buffer_end_index};
}

sub seq_buf_overrun
{
    my $self = shift;
    if (@_) { $self->{seq_buf_overrun} = shift; }
    return $self->{seq_buf_overrun};
}

sub seekpos {
    my $self = shift;
    my $pos = shift;
    
    seek($self->{FILE_H}, $pos, 0);
}

sub open_file
{
    my $self = shift;
    my $file = shift;
    my $mode = shift;
    
    my $success = 0;
    
    if ($mode eq "read") {
        &open_for_read(\$self->{FILE_H}, $file);
        $self->{seq_id} = 0;
        $self->{saved_line} = "";
    }
    elsif ($mode eq "write") {
        &open_for_write(\$self->{FILE_H}, $file);        
    }
    elsif ($mode eq "append") {
        &open_for_append(\$self->{FILE_H}, $file);        
    }
    $self->{file_name} = $file;
    $success = 1;

    return $success;
}

sub close_file
{
    my $self = shift;
    
    if (defined $self->{FILE_H}) {
        close($self->{FILE_H});
    }
}

# Reads length of sequence first, then pre-extends to total length
# before reading it in (important optimization for very long sequences)
# Also, will search for sequence name matching $key

sub read_fasta {
    
    my $self = shift;
    my $opts = shift;
    my $target_seq_id = shift;
    my $key = $opts->seq_key();
    my $fh = $self->{FILE_H};
    
    my ($seqlen, $filepos, $pre_extend_len, $seq_index_step, @seq_index);

    $self->{seq_name} = "";
    $self->{seq_description} = "";
    $self->{seq_length} = 0;
    $self->{sequence} = "";

# if $key is not the global $seq_key (non-alphanumerics already
#  escaped out for $seq_key) then escape out '\' problem causing char's
#    if ($key ne $seq_key) {
#        $key =~ s/(\W)/\\$1/g;
#    }        
    
    while ((!eof($fh)) 
           && (($self->{saved_line} =~ /^>/) || ($self->{saved_line} = <$fh>))) 
    {                                
        if (($self->{saved_line} =~ /^>\s*($key)\s+(.*)$/) ||
            ($opts->start_at_key()) && ($self->{key_found}) &&
            ($self->{saved_line} =~ /^>\s*(\S*)\s+(.*)$/o))
        {
            $self->{seq_id}++;

            # if searching for a particular SeqID go on to next seq
            #  if target and current seqid's don't match
            if ($target_seq_id && ($self->{seq_id} != $target_seq_id)) {
                $self->{saved_line} = <$fh>;
                next;
            }

            $self->{key_found}       = 1;
            $self->{seq_name}        = $1;
            $self->{seq_description} = $2;
            $self->{sequence}        = "";
            $self->{seq_name_map}->{$self->{seq_name}} = $self->{seq_id};
            
            @seq_index        = ();
            $seq_index_step   = $self->{seq_index_inc};   # set first bp position to save

            $filepos = tell($fh);
            $seqlen = 0;
            push(@seq_index, $seqlen, tell($fh));
            $pre_extend_len = 0;
#            print LOGFILE "At pos: ";

            while ($self->{saved_line} = <$fh>)
            {
                if ($self->{saved_line} =~ /^>/) { last; }
                $self->{saved_line} =~ s/[ \n\t\d]//g;     # strip whitespace & numbers
                $seqlen += length($self->{saved_line});
                
                # Save the start position of this chunk of seq for later easy return
                if ($seqlen > $seq_index_step) {
                    push(@seq_index, $seqlen, tell($fh));
                    $seq_index_step += $self->{seq_index_inc};
#                    print LOGFILE "($Seqlen) ";
                } 
                
                if (($pre_extend_len == 0) && ($seqlen >= $self->{max_seq_buffer})) {
                    $pre_extend_len = $seqlen;
                }
            }
            push(@seq_index, $seqlen, tell($fh));                        
            $self->{seq_length} = $seqlen;
#            print LOGFILE " ";
            
            $self->{all_seq_indices}->[$self->{seq_id}] = [@seq_index];

            seek($fh,$filepos,0);
            $self->{sequence} = 'X' x $pre_extend_len;  # pre-extending string for efficiency
            $seqlen = 0;
            while (($seqlen < $self->{max_seq_buffer}) && ($self->{saved_line} = <$fh>))
            {
                if ($self->{saved_line} =~ /^>/) { last; }
                $self->{saved_line} =~ s/[ \n\t\d]//g;     # strip whitespace & numbers
                substr($self->{sequence}, $seqlen, length($self->{saved_line})) = $self->{saved_line};
                $seqlen += length($self->{saved_line});
            }                        

            # if sequence is longer than MaxSeqBuffer length,
            # then save last ~200 nt to allow overlap with next buffer frame 
            # this prevents tRNAs on the border between buffers from being chopped
            # in half (and missed!)

            if ($seqlen >= $self->{max_seq_buffer}) {
                $self->{buffer_overlap_seq} = substr($self->{sequence}, $seqlen - $self->{seq_buf_overlap});
                $self->{buffer_end_index}   = $seqlen - length($self->{buffer_overlap_seq});
                $self->{seq_buf_overrun} = 1;
            }
            else {
                $self->{seq_buf_overrun} = 0;
            }
            
            $self->{buffer_length} = length($self->{sequence});
            $self->{sequence} = uc($self->{sequence});
            $self->{sequence} =~ s/U/T/g;
            $self->{sequence} =~ s/X/N/g;
            
            ## Remove long runs of N's from consideration by pre-scanners
            ## By doing this, pre-scanner false-pos rate is normal, even
            ## when scanning unfinished genomes with long N insert "placeholders"
            $self->{sequence} =~ s/NNNNNNNNNN/CCCCCCCCCC/g; 

            return 1;
        }
        else {
            if ($self->{saved_line} =~ /^>/) {
                $self->{seq_id}++;
            }
            $self->{saved_line} = <$fh>;
        }
    }                                
    0;                                
}
                
sub read_fasta_subseq  {
    
    my $self = shift;
    my $target_seq_id = shift;
    my $subseq_start = shift;
    my $subseq_len = shift;
    my $fh = $self->{FILE_H};
    
    my ($seqlen, $filepos, $curpos, $tempseq, $seq_head, $index_pos, $ct);

    $self->{seq_length} = 0;
    $self->{sequence} = "";

    # find closest position in desired sequence from file position index

    $ct=0;
    if (!defined $self->{all_seq_indices}->[$target_seq_id]) {
        $seqlen = 0;
        $index_pos = 0;
    }
    else {
        while ($self->{all_seq_indices}->[$target_seq_id][$ct] < $subseq_start) {
            $ct+=2;
        }
        $seqlen     = $self->{all_seq_indices}->[$target_seq_id][$ct-2]; 
        $index_pos  = $self->{all_seq_indices}->[$target_seq_id][$ct-1];
    }
    seek ($fh, $index_pos, 0);

    $tempseq = "";

    # scan until I get to the sequence position 

    while (($seqlen < $subseq_start) && ($self->{saved_line} = <$fh>))
    {
        if ($self->{saved_line} =~ /^>/) { 
            return 0; 
        }
        $self->{saved_line} =~ s/[ \n\t\d]//g;     # strip whitespace & numbers
        $seqlen += length($self->{saved_line});
    }

    $tempseq = 'X' x $subseq_len;  # pre-extending string for efficiency
            
    $curpos = $seqlen - length($self->{saved_line});
    $seq_head = substr($self->{saved_line}, $subseq_start - $curpos - 1); 
    substr($tempseq, 0, length($seq_head)) = $seq_head;
        
    $seqlen = length($seq_head);
            
    while (($seqlen < $subseq_len) && ($self->{saved_line} = <$fh>))
    {
        if ($self->{saved_line} =~ /^>/) { last; }
        $self->{saved_line} =~ s/[ \n\t\d]//g;     # strip whitespace & numbers
        substr($tempseq, $seqlen, length($self->{saved_line})) = $self->{saved_line};
        $seqlen += length($self->{saved_line});
    }                        
    
    $self->{sequence} = substr($tempseq, 0, $subseq_len);

    $self->{sequence} = uc($self->{sequence});
    $self->{sequence} =~ s/U/T/g;
    $self->{sequence} =~ s/X/N/g;
    $self->{seq_length} = length($self->{sequence});
    return 1;
}

sub read_fasta_subseq_slow {
    
    my $self = shift;
    my $opts = shift;
    my $key = shift;
    my $target_seq_id = shift;
    my $subseq_start = shift;
    my $subseq_len = shift;
    my $fh = $self->{FILE_H};
    
    my ($seqlen, $filepos, $curpos, $tempseq);
    my $last_header = "";
    my $seq_head = "";

    $self->{seq_length} = 0;
    $self->{sequence} = "";

# if $key is not the global $seq_key (non-alphanumerics already
#  escaped out for $seq_key) then escape out '\' problem causing char's
#  if ($key ne $seq_key) {
        $key =~ s/(\W)/\\$1/g;
#    }        

    while ((!eof(FAHANDLE)) 
           && (($self->{saved_line} =~ /^>/) || ($self->{saved_line} = <FAHANDLE>))) 
    {                                
        if (($self->{saved_line} =~ /^>\s*($key)\s+(.*)$/) ||
            ($opts->start_at_key()) && ($self->{key_found}) &&
            ($self->{saved_line} =~ /^>\s*(\S*)\s+(.*)$/o))
        {
            $self->{seq_id}++;
            
            # if searching for a particular SeqID go on to next seq
            #  if target and current seqid's don't match
            if ($target_seq_id && ($self->{seq_id} != $target_seq_id)) {
                $self->{saved_line} = <$fh>;
                next;
            }

            $filepos = tell($fh);  # save position of last fasta header
            $last_header = $self->{saved_line}; 
            
            $self->{key_found} = 1;
            $self->{seq_name}        = $1;
            $self->{seq_description} = $2;
            $self->{sequence}        = "";
            $tempseq = "";

            $seqlen = 0;
            while (($seqlen < $subseq_start) && ($self->{saved_line} = <$fh>)) {
                if ($self->{saved_line} =~ /^>/) { last; }
                $self->{saved_line} =~ s/[ \n\t\d]//g;     # strip whitespace & numbers
                $seqlen += length($self->{saved_line});
            }

            $tempseq = 'X' x $subseq_len;  # pre-extending string for efficiency
            
            $curpos = $seqlen - length($self->{saved_line});
            $seq_head = substr($self->{saved_line}, $subseq_start - $curpos - 1); 
            substr($tempseq, 0, length($seq_head)) = $seq_head;
        
            $seqlen = length($seq_head);
            
            while (($seqlen < $subseq_len) && ($self->{saved_line} = <$fh>)) {
                if ($self->{saved_line} =~ /^>/) { last; }
                $self->{saved_line} =~ s/[ \n\t\d]//g;      # strip whitespace & numbers
                substr($tempseq, $seqlen, length($self->{saved_line})) = $self->{saved_line};
                $seqlen += length($self->{saved_line});
            }                        

            $self->{sequence} = substr($tempseq, 0, $subseq_len);

            $self->{sequence} = uc($self->{sequence});
            $self->{sequence} =~ s/U/T/g;
            $self->{sequence} =~ s/X/N/g;
            $self->{seq_length} = length($self->{sequence});
            seek($fh, $filepos, 0);                 # return file position to beginning of this seq
            $self->{seq_id}--;                      # rewind seqid by 1
            $self->{saved_line} = $last_header;             # restore to original seq header line
            return 1;
        }
        else {
            if ($self->{saved_line} =~ /^>/) {
                $self->{seq_id}++;
            }
            $self->{saved_line} = <$fh>;
        }
    }                                
    0;                                
}

## read_more_fasta  
## Reads remaining portion of large fasta file (size>$MaxSeqBuffer)
## Only reads in $MaxSeqBuffer amount or less each time
                
sub read_more_fasta {
    
    my $self = shift;
    my $fh = $self->{FILE_H};
    
    my ($seqlen, $filepos);
    
    $filepos = tell($fh);
    $seqlen = 0;
    while (($seqlen + $self->{seq_buf_overlap} < $self->{max_seq_buffer}) && ($self->{saved_line} = <$fh>))
    {
        if ($self->{saved_line} =~ /^>/) { last; }
        $self->{saved_line} =~ s/[ \n\t\d]//g;     # strip whitespace & numbers
        $seqlen += length($self->{saved_line});
    }                        

    if ($seqlen == 0) {
        return 0;
    }

    seek($fh, $filepos, 0);

    $self->{sequence} = $self->{buffer_overlap_seq}. 'X' x $seqlen;  # pre-extending string for efficiency
    $seqlen = length($self->{buffer_overlap_seq});    

    while (($seqlen < $self->{max_seq_buffer}) && ($self->{saved_line} = <$fh>))
    {
        if ($self->{saved_line} =~ /^>/) { last; }
        $self->{saved_line} =~ s/[ \n\t\d]//g;     # strip whitespace & numbers
        substr($self->{sequence}, $seqlen, length($self->{saved_line})) = $self->{saved_line};
        $seqlen += length($self->{saved_line});
    }                        
    
    # if sequence is longer than MaxSeqBuffer length,
    # then save last ~200 nt to allow overlap with next buffer frame 
    # this prevents tRNAs on the border between buffers from being chopped
    # in half (and missed!)
    
    if ($seqlen >= $self->{max_seq_buffer}) {
        $self->{buffer_overlap_seq} = substr($self->{sequence}, $seqlen - $self->{seq_buf_overlap});
        $self->{buffer_end_index}   += $seqlen - length($self->{buffer_overlap_seq});
        $self->{seq_buf_overrun} = 1;
    }
    else {
        $self->{seq_buf_overrun} = 0;
    }
    
    $self->{buffer_length} = length($self->{sequence});
    $self->{sequence} = uc($self->{sequence});
    $self->{sequence} =~ s/U/T/g;
    $self->{sequence} =~ s/X/N/g;
    
    ## Remove long runs of N's from consideration by pre-scanners
    ## By doing this, pre-scanner false-pos rate is normal, even
    ## when scanning unfinished genomes with long N insert "placeholders"
    $self->{sequence} =~ s/NNNNNNNNNN/CCCCCCCCCC/g; 
    
    return 1;
}
        
sub write_fasta {
    
    my $self = shift;
    my $fh = $self->{FILE_H};
    
    my ($pos, $line);

    print $fh ">$self->{seq_name} $self->{seq_description}\n"; 
    for ($pos = 0; $pos < length($self->{sequence}); $pos += 60)
    {
        $line = substr($self->{sequence}, $pos, 60);
        print $fh $line, "\n";
    }
}

sub get_tRNA_sequence {
    
    my $self = shift;
    my ($src_seq_name, $strand, $start, $end, $log, $opts, $constants) = @_;
    
    $self->{seq_name} = $src_seq_name;
    $self->{seq_description} = ""; 
    
    my ($upstream_len, $downstream_len, $src_seq_len, $fwd_start, $query_len, $upstream, $downstream, $tRNA_seq);
    my $src_seqid = $self->get_seq_id_from_name($src_seq_name);
    
    $upstream_len   = $constants->upstream_len();
    $downstream_len = $constants->downstream_len();
    if ($strand) {
        if ($start - $upstream_len <= 0) {
            $upstream_len = $start - 1;
        }
        $fwd_start = $start - $upstream_len;
        $src_seq_len = $end - $start + 1;
    }
    else {
        if ($end - $downstream_len <= 0) {
            $downstream_len = $end - 1;
        }
        $fwd_start = $end - $downstream_len;
        $src_seq_len = $start - $end + 1;
    }
    $query_len = $upstream_len + $src_seq_len + $downstream_len;

    if (!$self->read_fasta_subseq($src_seqid, $fwd_start, $query_len)) {
        
        # if can't find it on first try, reposition 
        # to beginning of file & try once more
            
        $log->write_line("Missed $src_seq_name using quick index. Rewinding seq file and trying again with slow search...");
        $self->seekpos(0);
        if (!$self->read_fasta_subseq_slow($opts, $src_seq_name, $src_seqid, $fwd_start, $query_len)) {
            print STDERR "Could not find $src_seq_name in ".$opts->fastafile()."\n";
            $log->write_line("Skipping to next tRNA hit...");
            return 0;
        }
    }
    
    if ($strand) {
        $downstream_len = $self->{seq_length} - $upstream_len - $src_seq_len;
        $upstream = substr($self->{sequence}, 0, $upstream_len);
        $downstream = "";
        if ($downstream_len > 0) {
            $downstream = substr($self->{sequence}, $upstream_len + $src_seq_len);
        }
        $tRNA_seq = substr($self->{sequence}, $upstream_len, $src_seq_len);
    }
    else {
        $upstream_len = $self->{seq_length} - $downstream_len - $src_seq_len;        
        $self->{sequence} = &rev_comp_seq($self->{sequence});
        $upstream = "";
        if ($upstream_len > 0) {
            $upstream = substr($self->{sequence}, 0, $upstream_len);
        }
        $downstream = substr($self->{sequence}, $upstream_len + $src_seq_len);
        $tRNA_seq = substr($self->{sequence}, $upstream_len, $src_seq_len);        
    }
    return ($tRNA_seq, $upstream, $downstream);
}

sub mask_out_sequence {

    my $self = shift;
    my ($seq_file, $temp_seq_file, $r_sorted_cms_hits) = @_;
    
    my $cms_hit = undef;
    my $fh_seq_in = undef;
    my $fh_seq_out = undef;
    my $line = "";
    my $last_line = "";
    my $seqname = "";
    my %cms_hits = ();
    my $hits = [];
    my $ct = 0;
    my $written_len = 0;
    my $seq_start = 0;
    my $subseq_start = 0;
    my $subseq_len = 0;
    my $N_start = 0;
    
    foreach $cms_hit (@$r_sorted_cms_hits) {
        if (defined $cms_hits{$cms_hit->{seqname}}) {
            push (@{$cms_hits{$cms_hit->{seqname}}}, $cms_hit);
        }
        else {
            $hits = [];
            push (@$hits, $cms_hit);
            $cms_hits{$cms_hit->{seqname}} = $hits;
        }
    }
    
    &open_for_read(\$fh_seq_in, $seq_file);
    &open_for_write(\$fh_seq_out, $temp_seq_file);
    
    while ($line = <$fh_seq_in>) {
        chomp($line);
        if ($line =~ /^>([^\t]+)$/) {
            $seqname = $1;
            $seqname = &trim($seqname);
            if (index($seqname, ' ') > -1) {
                $seqname = substr($seqname, 0, index($seqname, ' '));
            }
            $hits = undef;
            $ct = 0;
            if (defined $cms_hits{$seqname}) {
                $hits = $cms_hits{$seqname};
                $subseq_len = $hits->[$ct]->{len};
                $seq_start = $hits->[$ct]->{start};
                $seq_start = $hits->[$ct]->{end} if ($hits->[$ct]->{strand} == 0); 
            }
            $written_len = 0;
            $N_start = 0;
            if ($last_line ne "") {
                print $fh_seq_out $last_line . "\n";
                $last_line = "";
            }
            print $fh_seq_out $line . "\n";
        }
        elsif ($line =~ /^\s*$/) {
        }
        else {
            if ($last_line ne "") {
                $line = $last_line . $line;
                $last_line = "";
            }
            if (defined $hits) {
                if ($ct < scalar(@$hits)) {
                    if ($written_len + length($line) < $seq_start) {
                        print $fh_seq_out $line . "\n";
                        $written_len += length($line);
                    }
                    else {
                        if ($N_start > 0) {
                            $subseq_start = 1;
                        }
                        else {
                            $subseq_start = $seq_start - $written_len;
                            print $fh_seq_out substr($line, 0, $subseq_start - 1);
                            $written_len += ($subseq_start - 1);
                        }
                        if (length($line) >= ($subseq_start + $subseq_len - 1)) {
                            print $fh_seq_out 'N' x $subseq_len;
                            print $fh_seq_out "\n";
                            $written_len += $subseq_len;
                            if (length($line) > ($subseq_start + $subseq_len - 1)) {
                                $last_line = substr($line, $subseq_start + $subseq_len - 1);
                            }
                            $N_start = 0;
                            $ct++;
                            if ($ct < scalar(@$hits)) {
                                $subseq_len = $hits->[$ct]->{len};
                                $seq_start = $hits->[$ct]->{start};
                                $seq_start = $hits->[$ct]->{end} if ($hits->[$ct]->{strand} == 0);
                            }
                        }
                        else {
                            print $fh_seq_out 'N' x (length($line) - $subseq_start + 1);
                            print $fh_seq_out "\n";
                            $written_len += (length($line) - $subseq_start + 1);
                            $N_start = 1;
                            $subseq_len -= (length($line) - $subseq_start + 1);
                        }
                    }
                }
                else {
                    print $fh_seq_out $line . "\n";
                    $written_len += length($line);
                }
            }
            else {
                print $fh_seq_out $line . "\n";
                $written_len += length($line);
            }
        }
    }
    
    close($fh_seq_in);
    close($fh_seq_out);
}

1;
