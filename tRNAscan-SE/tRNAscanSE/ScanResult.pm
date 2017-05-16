# tRNAscanSE/ScanResult.pm
# This class describes the outputs of scan results used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::ScanResult;

use strict;
use tRNAscanSE::Utils;
use tRNAscanSE::Sequence;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(init_fp_result_file save_Acedb_from_firstpass save_firstpass_output 
				prep_for_secpass_only parse_tabular_output write_tRNA output_tRNA output_split_fragments);

our $printed_header = 0;            # keeps track of whether or
                                    # or not results column header
                                    # has been printed yet
our ($max_seq_name_width, $max_seq_len_width);

sub init_fp_result_file {
    
    my ($file) = @_;
    
    &open_for_append(\*FILE_H, $file);
    
    print FILE_H "Sequence\t\ttRNA Bounds\ttRNA\tAnti\t\n";
	print FILE_H "Name     \ttRNA #\tBegin\tEnd\tType\tCodon\t",
	    "SeqID\tSeqLen\tScore\n";
	print FILE_H "--------\t------\t-----\t---\t----\t-----\t",
	    "-----\t------\t-----\n";
    
    close(FILE_H);
}

sub save_Acedb_from_firstpass  {

    my ($output_codon, $r_one_let_trans_map, $r_hit_list, $out_file) = @_;
    my($i, $triplet);

    &open_for_append(\*FILE_H, $out_file);

    foreach $i (0..(scalar(@$r_hit_list) - 1)) {
	printf FILE_H "Sequence\t%s\nSubsequence\t%s.t%d %d %d\n\n",
		$r_hit_list->[$i]{seqname}, $r_hit_list->[$i]{seqname},
		$i + 1, $r_hit_list->[$i]{start}, $r_hit_list->[$i]{end};
	
	printf FILE_H "Sequence\t%s.t%d\nSource\t\t%s\n",
		$r_hit_list->[$i]{seqname}, $i + 1, $r_hit_list->[$i]{seqname};
	if ($r_hit_list->[$i]{istart} > 0) {
	    if ($r_hit_list->[$i]{istart} < $r_hit_list->[$i]{iend}) {
		printf FILE_H "Source_Exons\t1 %d\n",
			$r_hit_list->[$i]{istart} - $r_hit_list->[$i]{start};
		printf FILE_H "Source_Exons\t%d %d\n",
			$r_hit_list->[$i]{iend} - $r_hit_list->[$i]{start} + 2,
			$r_hit_list->[$i]{end} - $r_hit_list->[$i]{start} + 1; }
	    else {
		printf FILE_H "Source_Exons\t1 %d\n",
			$r_hit_list->[$i]{start} - $r_hit_list->[$i]{istart} + 1;
		printf FILE_H "Source_Exons\t%d %d\n",
			$r_hit_list->[$i]{start} - $r_hit_list->[$i]{iend} + 2,
			$r_hit_list->[$i]{start} - $r_hit_list->[$i]{end} + 1; }
	}	 
	printf FILE_H "Brief_identification tRNA-%s\n", $r_hit_list->[$i]{type};
	
	# either output Codon or Anticodon for tRNA
	$triplet = uc($r_hit_list->[$i]{acodon});
	if ($output_codon) {
	    $triplet = &rev_comp_seq($triplet);
	}

	printf FILE_H "Transcript tRNA \"%s %s %s\"\n\n",
	    $triplet, $r_hit_list->[$i]{type}, $r_one_let_trans_map->{$r_hit_list->[$i]{type}};
	
    }
    close(FILE_H);
}

sub print_results_header {
    
    my ($opts, $get_hmm_score, $seq_name_width, $seq_len_width) = @_;
    my ($label, $codon_label) = "";
    
    if ($opts->cove_mode()) {
		$label = "\tCove";
    }
    elsif ($opts->infernal_mode()) {
		$label = "\tCM";
    }
    elsif ($opts->eufind_mode() && !$opts->tscan_mode()) {
		$label = "\tEufind";
    }

    if ($opts->output_codon()) {
		$codon_label = "   "; 
    }
    else {
		$codon_label = "Anti";
    }
    
    if (!($opts->ace_output())) {
		&open_for_append(\*OUTFILE, $opts->out_file());
	
		printf OUTFILE "%-".$seq_name_width."s\t\t","Sequence";
		printf OUTFILE "%-".$seq_len_width."s\t","tRNA";
		printf OUTFILE "%-".$seq_len_width."s\t","Bounds";
		print  OUTFILE "tRNA\t$codon_label\tIntron Bounds",$label;
	
		if  ($get_hmm_score) { 
			print OUTFILE "\tHMM\t2'Str";
		}
		if ($opts->save_source()) {
			print OUTFILE "\tHit";
		}
		if ($opts->search_mode() eq "archaea") {
			print OUTFILE "\tIntron";
		}
		print OUTFILE "\n";

		printf OUTFILE "%-".$seq_name_width."s\t","Name";
		print  OUTFILE "tRNA \#\t";
		printf OUTFILE "%-".$seq_len_width."s\t","Begin";
		printf OUTFILE "%-".$seq_len_width."s\t","End";
	
		print OUTFILE "Type\tCodon\tBegin\tEnd\tScore";
	
		if  ($get_hmm_score) { 
			print OUTFILE "\tScore\tScore";
		}
		if ($opts->save_source()) {
			print OUTFILE "\tOrigin";
		}
		if ($opts->search_mode() eq "archaea") {
			print OUTFILE "\tCount";
		}
		print OUTFILE "\n";
	
		printf OUTFILE "%-".$seq_name_width."s\t","--------";
		print  OUTFILE "------\t";
		printf OUTFILE "%-".$seq_len_width."s\t","----";
		printf OUTFILE "%-".$seq_len_width."s\t","------";
		print  OUTFILE "----\t-----\t-----\t----\t------";
	
		if  ($get_hmm_score) { 
			print OUTFILE "\t-----\t-----";
		}
		if ($opts->save_source()) {
			print OUTFILE "\t------";
		}
		if ($opts->search_mode() eq "archaea") {
			print OUTFILE "\t----------";
		}
		print OUTFILE "\n";
    }
    close OUTFILE;
}

sub save_firstpass_output {

    my ($opts, $r_hit_list, $r_source_tab, $r_fpass_trna_base_ct, $seq_len, $seq_id) = @_;
    my ($i, $triplet);
    
    if (!$opts->CM_mode()) {
		if (!($opts->brief_output() || $printed_header)) {
			&print_results_header($opts, 0, 20, 20);
			$printed_header = 1;
		}
		&open_for_append(\*TAB_RESULTS, $opts->out_file());
    }
    else {		       
		&open_for_append(\*TAB_RESULTS, $opts->firstpass_result_file());	
    }
    
    foreach $i (0..(scalar(@$r_hit_list) - 1)) {

		$triplet = uc($r_hit_list->[$i]{acodon});
		if ($opts->output_codon()) {
			$triplet = &rev_comp_seq($triplet);
		}
		
		printf TAB_RESULTS "%-10s\t%d\t%d\t%d\t%s\t%s\t",
			$r_hit_list->[$i]{seqname}, $i + 1,
			$r_hit_list->[$i]{start}, $r_hit_list->[$i]{end},
			$r_hit_list->[$i]{type}, $triplet;
		
		# save intron bounds if not doing Cove analysis
		
		if (!$opts->CM_mode()) {
			printf TAB_RESULTS "%d\t%d\t%.2f", $r_hit_list->[$i]{istart},
			$r_hit_list->[$i]{iend}, $r_hit_list->[$i]{score};
		}
	
		# save seq id number and source seq length if needed for Cove analysis 
	
		else {
			printf TAB_RESULTS "%d\t%d\t%.2f", $seq_id, $seq_len, $r_hit_list->[$i]{score};
		}
		
		if ($opts->save_source()) {
			print TAB_RESULTS " ", $r_source_tab->[$r_hit_list->[$i]{source}];
		}
		print TAB_RESULTS "\n";
		
		$$r_fpass_trna_base_ct += abs($r_hit_list->[$i]{end} - $r_hit_list->[$i]{start}) + 1;
    }
    close TAB_RESULTS;
}				

# Create dummy first-pass result file with all sequences
sub prep_for_secpass_only  {       

    my ($opts, $stats, $seq_file) = @_;
    my ($saved_line, $targ_seq_id);

    $seq_file->open_file($opts->fasta_file(), "read");

    &open_for_append(\*RESFILE, $opts->firstpass_result_file());	
    $saved_line = '';
    $targ_seq_id = 0;      # Don't look for a specific Seq number
 
    while ($seq_file->read_fasta($opts, $targ_seq_id)) {
		print (RESFILE $seq_file->seq_name()."\t1\t1\t".$seq_file->seq_length()."\t???\t???\t".$seq_file->seq_id()."\t".$seq_file->seq_length()." C\n");
		print (RESFILE $seq_file->seq_name()."\t2\t".$seq_file->seq_length()."\t1\t???\t???\t".$seq_file->seq_id()."\t".$seq_file->seq_length()." C\n");

		$stats->increment_numscanned();
    }
    close RESFILE;
    $seq_file->close_file();
}

# read first pass result file one input sequence at a time,
# putting results in array @prescan_tRNAs

sub parse_tabular_output {

    my ($opts, $r_prescan_tRNAs, $r_seqinfo_flag) = @_;
    my $firstpass_result_file = $opts->firstpass_result_file();
    my $padding = $opts->padding();
    
    my ($seq_name, $trnact, $trnaName,
        $ts_start, $ts_end, $ts_len, $sense_strand,
        $ts_seq_id, $ts_seq_len, $score, $ts_type, $ts_anticodon,
	$hit_source);
    
    # open first-pass tabular result file
    open (FIRSTPASS_TRNAS, "$firstpass_result_file") || 
	die "FATAL: Can't open first-pass tRNA output file $firstpass_result_file\n\n" ; 
    
    while (<FIRSTPASS_TRNAS>) 
    {	
		if (/Type\tCodon\tSeqID\tSeqLen/)  {
			# Column header present if we record seqID's and lengths
			$$r_seqinfo_flag = 1;
		}
		elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/o)  
		{
			$seq_name = $1;
			$trnact  = $2;
			$trnaName = $seq_name.".t".$trnact;
			$ts_start = $3;	        # trna subseq absolute start index
			$ts_end = $4;		# trna subseq absolute end index
			$ts_type = $5;
			$ts_anticodon = $6;
			$ts_seq_id = $7;
			$ts_seq_len = $8;
			$score = $9;
			$hit_source = $';
			$hit_source =~ s/[\s\t\n]//g; 
			
			# if seqinfo_flag not set, file does not have SeqID info in
			#  7th column of output, don't mistake number read for SeqID
			
			if (!$$r_seqinfo_flag) {
				$ts_seq_id = 0;
			}
			
			if ($ts_end > $ts_start)  
			{
				$sense_strand = 1;     # flag for forward or reverse strand
			
				# pad ends of sequence only if EufindtRNA is being used
				#  and $seqinfo_flag is set (we know the seq lengths)
	
				if ($opts->eufind_mode() && $$r_seqinfo_flag) 
				{
					$ts_start = &max(1, $ts_start - $padding);
					$ts_end   = &min($ts_seq_len, $ts_end + $padding);
				}
				$ts_len = $ts_end - $ts_start + 1;
			}
			else  { 
				$sense_strand = 0;
				if ($opts->eufind_mode() && $$r_seqinfo_flag) {
					$ts_start = &min($ts_seq_len, $ts_start + $padding);
					$ts_end = &max(1, $ts_end - $padding);
				}
				$ts_len = $ts_start - $ts_end + 1;
		    }
	    
			if ($ts_end == $ts_start) {
				print STDERR "Error reading $firstpass_result_file: tRNA of length 0"; 
			}
	    
			push(@$r_prescan_tRNAs,
			 {seq => "", name => $trnaName, 
			  start => $ts_start, end => $ts_end, len => $ts_len, 
			  isotype => $ts_type, acodon => $ts_anticodon, score => $score,
			  src_seqname => $seq_name, src_seqlen => $ts_seq_len, 
			  src_seqid => $ts_seq_id, strand => $sense_strand, 
			  hit_source => $hit_source});
			
		}	# while <FIRSTPASS_TRNAS> not at eof
    }
 
    close FIRSTPASS_TRNAS;
}

sub write_tRNA {
    
    my ($file_name, $seq_name, $seq_desc, $seq, $overwrite) = @_;
    
    my $trna_file = tRNAscanSE::Sequence->new;
    my $write_mode = "append";
    if ($overwrite) {
		$write_mode = "write";
    }
    $trna_file->set_seq_info($seq_name, $seq_desc, length($seq), $seq);
    $trna_file->open_file($file_name, $write_mode);
    $trna_file->write_fasta();
    $trna_file->close_file();
}

# Write final tRNA prediction to various selected output sources/files
# Sets globals $MaxSeqNameWidth and $MaxSeqLenWidth and $printed_header

sub output_tRNA {

    my ($opts, $gc, $log, $r_tab_results, $get_hmm_score, $program_id,
	$r_fp_tRNA_info,           # first pass scanner tRNA info 
	$r_tRNA_info,              # final tRNA info
	$curseq_trnact) = @_;

    my $results_line = "";
    
    if (!$opts->results_to_stdout()) {
		$log->write_line("$r_tRNA_info->{ID}:  ".$opts->second_pass_label()." type= $r_tRNA_info->{isotype}\t ".
		"First-pass scan ($r_fp_tRNA_info->{hit_source}) type= $r_fp_tRNA_info->{isotype}\t".
		"Score= $r_tRNA_info->{score}");
    }
    if ($opts->save_all_struct()) {
		&save_allStruct_output($opts, $gc, $get_hmm_score, $r_tRNA_info, $curseq_trnact);
    }
    
    # Create tabular results line, ready for output
    
    if (!$printed_header) {
		$max_seq_name_width = &max(length($r_fp_tRNA_info->{src_seqname}) + 1, 8);
		$max_seq_len_width  = length($r_fp_tRNA_info->{src_seqlen});
    }
    
    $results_line = &construct_tab_output($opts, $get_hmm_score, $r_tRNA_info, $curseq_trnact, $max_seq_name_width, $max_seq_len_width);
    
    # Internal copy of results saved for later uses
    push(@$r_tab_results, $results_line);
    
    if ($opts->ace_output()) {       
		&save_Acedb_from_secpass($opts, $gc, $r_tRNA_info, $program_id);
    }
    else 
    {    
		if (!($opts->brief_output() || $printed_header)) {
			&print_results_header($opts, $get_hmm_score, $max_seq_name_width, $max_seq_len_width);
			$printed_header = 1;
	}	    
	&open_for_append(\*TABOUT, $opts->out_file());
	print TABOUT $results_line;
	close TABOUT;	    
	}			
}

sub save_allStruct_output {
    
    my ($opts, $gc, $get_hmm_score, $r_tRNA_info, $curseq_trnact) = @_;

    my $ruler = '    *    |' x 20;     # ruler printed out with
                                       #  secondary structure output

    my $seqlen = length($r_tRNA_info->{seq});

    &open_for_append(\*SECSTRUCT, $opts->all_struct_file());
    
    print SECSTRUCT "$r_tRNA_info->{seqname}.trna$curseq_trnact ($r_tRNA_info->{start}-$r_tRNA_info->{end})\t",
        "Length: $seqlen bp\nType: $r_tRNA_info->{isotype}\t";

    if ($opts->output_codon()) {
		print SECSTRUCT "Codon: ", &rev_comp_seq($r_tRNA_info->{acodon}), " at ";
    }
    else {
		print SECSTRUCT "Anticodon: $r_tRNA_info->{acodon} at ";
    }

    if ($r_tRNA_info->{acodon} eq $gc->undef_anticodon()) {
		print SECSTRUCT "0-0 (0-0)\t";
    }
    else {
		print SECSTRUCT "$r_tRNA_info->{acodon_pos}-", $r_tRNA_info->{acodon_pos} + 2;
		if (!$opts->arch_mode()) {
			if ($r_tRNA_info->{strand}) {
				print SECSTRUCT " (", $r_tRNA_info->{acodon_pos} + $r_tRNA_info->{start} - 1, "-",
						$r_tRNA_info->{acodon_pos} + $r_tRNA_info->{start} + 1,")\t";
			}
			else {
				print SECSTRUCT " (", $r_tRNA_info->{start} - $r_tRNA_info->{acodon_pos} + 1, "-",
						$r_tRNA_info->{start} - $r_tRNA_info->{acodon_pos} - 1,")\t";
			}
		}
		else {
				print SECSTRUCT " (", $r_tRNA_info->{acodon_pos}, "-", $r_tRNA_info->{acodon_pos} + 2,")\t";			
		}
	}

    print SECSTRUCT "Score: $r_tRNA_info->{score}\n";
    if (scalar(@{$r_tRNA_info->{introns}}) > 0) {
	
		foreach my $intron (@{$r_tRNA_info->{introns}}) {
			if (defined $intron) {
				if ($intron->{seq} ne "") {
					print SECSTRUCT "Possible intron: $intron->{start}-$intron->{end} ";
					if ($r_tRNA_info->{strand}) {	
						print SECSTRUCT "(", $intron->{start} + $r_tRNA_info->{start} - 1, "-",
							$intron->{end} + $r_tRNA_info->{start} - 1,")\n";
					}
					else {
						print SECSTRUCT "(", $r_tRNA_info->{start} - $intron->{start} + 1, "-",
							$r_tRNA_info->{start} - $intron->{end} + 1,")\n";
					}
				}
			}
		}
	}
	
    if ($r_tRNA_info->{is_pseudo}) {
		printf SECSTRUCT 
			"Possible pseudogene:  HMM Sc=%.2f\tSec struct Sc=%.2f\n",
			$r_tRNA_info->{hmm_score}, $r_tRNA_info->{ss_score};
    }
    elsif ($get_hmm_score) {
		printf SECSTRUCT 
			"HMM Sc=%.2f\tSec struct Sc=%.2f\n", $r_tRNA_info->{hmm_score}, $r_tRNA_info->{ss_score};
    }
    
    print SECSTRUCT "     ",substr($ruler, 0, $seqlen - 1),"\n";
    print SECSTRUCT "Seq: $r_tRNA_info->{seq}\nStr: $r_tRNA_info->{ss}\n";
	if (defined $r_tRNA_info->{precursor}) {
		foreach my $intron (@{$r_tRNA_info->{introns}}) {
			if (defined $intron) {
				my $intron_seq = uc($intron->{seq});
				if ($intron_seq ne "") {
					$r_tRNA_info->{precursor} =~ s/$intron_seq/\[$intron_seq\]/;
				}
			}
		}
		print SECSTRUCT "Pre: ". uc($r_tRNA_info->{precursor}) ."\n\n";		
	}
	else {
		print SECSTRUCT "\n";
	}
    close(SECSTRUCT);
}

# Save tRNA hits in Tabular output

sub construct_tab_output {

    my ($opts, $get_hmm_score, $r_tRNA_info, $curseq_trnact, $max_seq_name_width, $max_seq_len_width) = @_;
    
    my ($result_line, $tRNA_type);
    
    if ($r_tRNA_info->{is_pseudo}) {
		$tRNA_type = "Pseudo";
    }
    else {
		$tRNA_type = $r_tRNA_info->{isotype};
    }
    
    $result_line =  sprintf "%-".$max_seq_name_width."s\t", $r_tRNA_info->{seqname};
    $result_line .= "$curseq_trnact\t";
    
    $result_line .= sprintf "%-".$max_seq_len_width."d\t", $r_tRNA_info->{start};
    $result_line .= sprintf "%-".$max_seq_len_width."d\t", $r_tRNA_info->{end};
    
    $result_line .= "$tRNA_type\t";

    if ($opts->output_codon()) {
		$result_line .= (&rev_comp_seq($r_tRNA_info->{acodon}))."\t";
    }
    else {
		$result_line .= "$r_tRNA_info->{acodon}\t";
    }

    if (scalar(@{$r_tRNA_info->{introns}}) == 0) {
		$result_line .= "0\t0"; 
    }
    else {
		my $intron_ct = 0;
		for (my $i = 0; $i < scalar(@{$r_tRNA_info->{introns}}); $i++) {
			if (defined $r_tRNA_info->{introns}->[$i]) {
				if ($r_tRNA_info->{introns}->[$i]->{seq} ne "") {
					if ($intron_ct > 0) {
						$result_line .= ",";
					}
					if ($r_tRNA_info->{strand}) {	
						$result_line .= ($r_tRNA_info->{introns}->[$i]->{start} + $r_tRNA_info->{start}-1);
					}
					else {
						$result_line .= ($r_tRNA_info->{start} - $r_tRNA_info->{introns}->[$i]->{start}+1);
					}
					$intron_ct++;
				}
			}
		}
		$result_line .= "\t";
		$intron_ct = 0;
		for (my $i = 0; $i < scalar(@{$r_tRNA_info->{introns}}); $i++) {
			if (defined $r_tRNA_info->{introns}->[$i]) {
				if ($r_tRNA_info->{introns}->[$i]->{seq} ne "") {
					if ($intron_ct > 0) {
						$result_line .= ",";
					}
					if ($r_tRNA_info->{strand}) {	
						$result_line .= ($r_tRNA_info->{introns}->[$i]->{end} + $r_tRNA_info->{start} - 1);
					}
					else {
						$result_line .= ($r_tRNA_info->{start} - $r_tRNA_info->{introns}->[$i]->{end} + 1);
					}
					$intron_ct++;
				}
			}
		}
	}			
    $result_line .= "\t$r_tRNA_info->{score}";
 
    if ($get_hmm_score) {
		$result_line .= sprintf "\t%.2f\t%.2f", $r_tRNA_info->{hmm_score}, $r_tRNA_info->{ss_score};
    }
    if ($opts->save_source()) {
		$result_line .= "\t$r_tRNA_info->{hit_source}";
    }
	if ($opts->search_mode() eq "archaea") {
		if (scalar(@{$r_tRNA_info->{introns}}) == 0) {
			$result_line .= "\t"; 
		}
		else {
			my $ci_count = 0;
			my $nci_count = 0;
			for (my $i = 0; $i < scalar(@{$r_tRNA_info->{introns}}); $i++) {
				if ($r_tRNA_info->{introns}->[$i]->{type} eq "CI") {
					$ci_count++;
				}
				elsif ($r_tRNA_info->{introns}->[$i]->{type} eq "NCI") {
					$nci_count++;
				}
			}
			$result_line .= "\t";
			if ($ci_count > 0) {
				$result_line .= $ci_count . " CI";
			}
			if (($ci_count > 0) && ($nci_count > 0)) {
				$result_line .= " ";
			}
			if ($nci_count > 0) {
				$result_line .= $nci_count . " NCI";
			}
		}
	}
    $result_line .= "\n";
    
    return $result_line;
}

sub Save_Acedb_from_secpass {

    my ($opts, $gc, $r_tRNA_info, $program_id) = @_;                     

    &open_for_append(\*ACEOUT, $opts->out_file());

    print ACEOUT "Sequence\t$r_tRNA_info->{seqname}\nSubsequence\t$r_tRNA_info->{ID} $r_tRNA_info->{start} $r_tRNA_info->{end}\n\n";
    print ACEOUT "Sequence\t$r_tRNA_info->{ID}\nSource\t\t$r_tRNA_info->{seqname}\n";
    if ($r_tRNA_info->{iseq}) {
	print ACEOUT "Source_Exons\t1 ", $r_tRNA_info->{istart} - 1,"\n";
	print ACEOUT "Source_Exons\t", $r_tRNA_info->{iend} + 1," ", abs($r_tRNA_info->{end} - $r_tRNA_info->{start}) + 1,"\n";
    }	   
    print ACEOUT "Brief_identification tRNA-$r_tRNA_info->{isotype}\n",
        "Transcript tRNA \"";

    if ($opts->output_codon()) {
	print ACEOUT &rev_comp_seq($r_tRNA_info->{acodon});
    }
    else {
	print ACEOUT $r_tRNA_info->{acodon};
    }
    
    print ACEOUT " $r_tRNA_info->{isotype} ", $gc->one_let_trans_map()->{$r_tRNA_info->{isotype}},
        "\"\nScore $program_id $r_tRNA_info->{score}\n";

    if ($r_tRNA_info->{is_pseudo}) {
	printf ACEOUT "Remark \"Likely pseudogene (HMM Sc=%.2f / Sec struct Sc=%.2f)\"\n",
            $r_tRNA_info->{hmm_score},$r_tRNA_info->{ss_score};
    }
    print ACEOUT "\n";
    close ACEOUT;
}

sub output_split_fragments {

    my ($opts, $r_pairs, $r_5half_hits, $r_3half_hits) = @_;                     
	
	my ($r_5half, $r_3half);
	
	&open_for_append(\*SPLITFILE, $opts->split_fragment_file());
	printf SPLITFILE "Fragment1\tFragment2\tSeqName1\tStartPos1\tEndPos1\tSeqName2\tStartPos2\tEndPos2\tScore1\tScore2\n";
	
	foreach my $r_pair (@$r_pairs) {
		if (defined $r_pair->{"5h"} && defined $r_pair->{"3h"}) {
			$r_5half = $r_5half_hits->[$r_pair->{"5h"}];
			$r_3half = $r_3half_hits->[$r_pair->{"3h"}];
			print SPLITFILE $r_5half->{seq}."\t".$r_3half->{seq}."\t",
				$r_5half->{hit_seqname}."\t".$r_5half->{tRNA_start}."\t".$r_5half->{tRNA_end}."\t",
				$r_3half->{hit_seqname}."\t".$r_3half->{tRNA_start}."\t".$r_3half->{tRNA_end}."\t",
				$r_5half->{score}."\t".$r_3half->{score}."\n";
			print SPLITFILE $r_5half->{ss}."\t".$r_3half->{ss}."\t\t\t\t\t\t\t\t\n";
		}
		elsif (defined $r_pair->{"5h"} && !defined $r_pair->{"3h"}) {
			$r_5half = $r_5half_hits->[$r_pair->{"5h"}];
			print SPLITFILE $r_5half->{seq}."\t\t",
				$r_5half->{hit_seqname}."\t".$r_5half->{tRNA_start}."\t".$r_5half->{tRNA_end}."\t",
				"\t\t\t",
				$r_5half->{score}."\t\n";
			print SPLITFILE $r_5half->{ss}."\t\t\t\t\t\t\t\t\t\n";
		}
		elsif (!defined $r_pair->{"5h"} && defined $r_pair->{"3h"}) {
			$r_3half = $r_3half_hits->[$r_pair->{"3h"}];
			print SPLITFILE "\t".$r_3half->{seq}."\t",
				"\t\t\t",
				$r_3half->{hit_seqname}."\t".$r_3half->{tRNA_start}."\t".$r_3half->{tRNA_end}."\t",
				"\t".$r_3half->{score}."\n";
			print SPLITFILE "\t".$r_3half->{ss}."\t\t\t\t\t\t\t\t\n";
		}
	}
}

1;
