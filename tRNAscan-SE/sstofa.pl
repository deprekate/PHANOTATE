#! /usr/bin/perl


if ($#ARGV < 0) {
    die "\nConvert secondary structure file to fasta format\n",
	"\nUsage: sstofa <sec struct file>\n\n";
}				

$ss_file = shift;

open (SSFILE, $ss_file) || die "Couldn't find $ss_file\n";
open(SEQFILE,">-");
$ct = 0;

while ($line = <SSFILE>) {
    if ($line =~ /^(\S+)\s+(\(\d+\-\d+\))\s+Length:\s(\d+)\sbp/) {
	$SeqName = $1;
	$bounds = $2;
	$SeqLen = $3;
    }
    elsif ($line =~ /^Type:\s(\S+)\s+Anticodon:\s(\S+).+Score:\s(\S+)/) {
	$isotype = $1;
	$ac = $2;
	$SeqName .= "-".$isotype.$ac;
	$score = $3;
	$SeqDescription = "$bounds  $isotype ($ac) $SeqLen bp  Sc: $score";
    }
    elsif ($line =~ /pseudogene/) {
	$SeqDescription .= "  Pseudo";
    }
    elsif ($line =~ /^Seq:\s(\S+)$/) {
	$Seq = $1;
	&write_fasta($SeqName,$SeqDescription,length($Seq),
		     uc($Seq),SEQFILE);
	$SeqName = "";
	$bounds = "";
	$SeqLen = 0;
	$isotype = "";
	$ac = "";
	$score = 0.0;
	$SeqDescription = "";
	$Seq = "";
    }
}



# End Main   
    
sub write_fasta {
    local($name, $description, $length, $sequence,*FAHANDLE) = @_;
    local($pos, $line);

    print FAHANDLE ">$name $description\n"; 
    for ($pos = 0; $pos < $length; $pos += 60)
    {
	$line = substr($sequence,$pos,60);
	print FAHANDLE $line, "\n";
    }
    1;
}
	
