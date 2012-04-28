#!perl
#
# Description: Calculate allele frequency and genotype counts from testvariants data file.  Combine with gene annotation info (if available) into one file.
#
#
#
# Created by Jessica on 2012-04-09

use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);


my ($testvarfile, $outputfile, $geneannotationlistvar, %geneannotations, $dbnsfpfile);
my @dbnsfp_keepcol = c(14, 17, 21..34, 37);


GetOptions(
	'testvar=s' => \$testvarfile, 
	'o=s' => \$outputfile,
	'genelistvar=s' => \$geneannotationlistvar,
	'dbnsfp=s' => \$dbnsfpfile,
);

if (!defined $testvarfile) {
	optionUsage("option --testvar not defined\n");
} elsif (!defined $geneannotationlistvar) {
	optionUsage("option --genelistvar not defined\n");
}	elsif (!defined $outputfile) {
	optionUsage("option --o not defined\n");
} elsif (!defined $dbnsfpfile) {
	optionUsage("option --dbnsfp not defined\n");
}

# load in gene annotations from files
print STDERR "Preloading gene annotations from $geneannotationlistvar\n";
if ($geneannotationlistvar =~ /.bz2$/) {
	open (GENELISTVAR, "bzcat $geneannotationlistvar |") or die "Cannot read $geneannotationlistvar: $!\n";
} else {
	open (GENELISTVAR, "$geneannotationlistvar") or die "Cannot read $geneannotationlistvar: $!\n";
}
<GENELISTVAR>;
while (<GENELISTVAR>) {
	$_ =~ s/\s+$//;
	loadGeneAnnot($_, \%geneannotations);
}
close GENELISTVAR;


# prepare dbNSFP file for reading
open( my $dbnsfp_handle, "zcat $dbnsfpfile |" ) or die "Couldn't read $dbnsfpfile: $!";
my $dbnsfp_headerline = <$dbnsfp_handle>;
my @dbnsfp_header = split("\t", $dbnsfp_headerline);
my $dbnsfp_nextchrline;			# read through dbNSFP file until new chr, but that line needs to be saved
my %dbnsfp_annots;
my $dbnsfp_annot_ref = \%dbnsfp_annots;


open (OUT, ">$outputfile");
print OUT "chr\tbegin\tend\tvartype\tref\talt\txref\taltfreq\tMAF\tnhomref\tnhet\tnhomalt\tnmiss\tnhommaj\tnhommin\tCR\tgenesymbol\torientation\tcomponent\tcomponentIndex\thasCodingRegion\timpact\tnucleotidePos\tproteinPos\tannotationRefSequence\tsampleSequence\tgenomeRefSequence";
print OUT join("\t", @dbnsfp_header[@dbnsfp_keepcol])."\n";
if ($testvarfile =~ /.bz2$/) {
	open (FILE, "bzcat $testvarfile |") or die "Cannot read $testvarfile: $!\n";
} else {
	open (FILE, "$testvarfile") or die "Cannot read $testvarfile: $!\n";
}
open (FILE, "$testvarfile") or die "Cannot read $testvarfile file.\n";
my $header = <FILE>;
$header =~ s/\s+$//;
my @headerline = split("\t", $header);
my @ASMids = @headerline[8..$#headerline];
my $currchr = 'NA';
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($thischr, $thisstart, $thisend) = @line[1..3];
	my $vartype = $line[4];
	my ($ref, $alleleSeq) = @line[5..6];
	
	if ($thischr eq 'chr3') {
		exit;
	}
	
	if ($thischr ne $currchr) {
		print STDERR "Working on $thischr\n";			
		$currchr = $thischr;
		%dbnsfp_annots = ();							# clear out previous chromosome's annotations
		if ($currchr ne 'NA') {
			my @nextline = split("\t", $dbnsfp_nextchrline);
			$dbnsfp_annots{$nextline[1]} = @nextline[@dbnsfp_keepcol];		# load the first line of the new chr
		}
		($dbnsfp_nextchrline, $dbnsfp_annot_ref) = loaddbNSFP($dbnsfp_annot_ref, $thischr, $dbnsfp_handle, \@dbnsfp_keepcol);
	}

	my $thiscoord = join('_', @line[1..6]);
	
	my @genotypes = @line[8..$#line];
	my @genocounts = (0,0,0,0);
	foreach my $geno (@genotypes) {
		if ($thischr =~ 'chrX') {
			if ($geno eq '00')  { $genocounts[0]++; }
			if ($geno eq '01' || $geno eq '10')  { $genocounts[1]++; }
			if ($geno eq '11')  { $genocounts[2]++; }			
			if ($geno =~ 'N')  { $genocounts[3]++; }
		} elsif ($thischr =~ 'chrM' || $thischr =~ 'chrY') {
			if ($geno eq '0')  { $genocounts[0]++; }
			if ($geno eq '1' || $geno eq '10')  { $genocounts[1]++; }
			if ($geno =~ 'N')  { $genocounts[3]++; }
		} else {
			if ($geno eq '00')  { $genocounts[0]++; }
			if ($geno eq '01' || $geno eq '10')  { $genocounts[1]++; }
			if ($geno eq '11')  { $genocounts[2]++; }
			if ($geno =~ 'N')  { $genocounts[3]++; }
		}
	}
	
	my $ngenosubj = $genocounts[0]+$genocounts[1]+$genocounts[2];	
	
	if ($ngenosubj != 0) {	
		my $altfreq = ($genocounts[1]+2*$genocounts[2])/(2*$ngenosubj);
		my $callrate = 1-($genocounts[3]/scalar(@genotypes));
	
		my ($maf, $genocounts_relhutt_maj, $genocounts_relhutt_min) = (0,0,0);		# count number of each genotype class using the hutterite minor allele as the alternative allele
		if ($genocounts[0] >= $genocounts[2]) {			# more or same ref alleles than alt alleles
			$maf = ($genocounts[1]+2*$genocounts[2])/(2*$ngenosubj);
			$genocounts_relhutt_maj = $genocounts[0];
			$genocounts_relhutt_min = $genocounts[2];
		} else {																	# same or more alt alleles as ref alleles
			$maf = ($genocounts[1]+2*$genocounts[0])/(2*$ngenosubj);
			$genocounts_relhutt_maj = $genocounts[2];
			$genocounts_relhutt_min = $genocounts[0];
		}

		my $geneannot = (("\t") x 11);
		if (exists $geneannotations{$thiscoord}) {
			$geneannot = $geneannotations{$thiscoord};
		} 
		
		for (my $i=1; $i<=7; $i++) {
			if ($line[$i]) {
				print OUT "$line[$i]\t";
			} else {
				print OUT "\t";
			}
		}
		print OUT "$altfreq\t$maf\t".join("\t", @genocounts)."\t$genocounts_relhutt_maj\t$genocounts_relhutt_min\t$callrate\t$geneannot";
		
		my $dbnsfp_lookup = join("_", ($thischr, $thisend, $ref, $alleleSeq));
		if (($thisend-$thisstart)==1 && exists $dbnsfp_annots{$thisend}) {		# haven't decided what to print out for ins/del/complex variants that span multiple bases
			print OUT "$dbnsfp_annots{$thisend}\n";		
		} else {
			print OUT (("\t") x scalar(@dbnsfp_keepcol))."\n";			
		}
	}
}
close FILE;
close OUT;

close $dbnsfp_handle;


sub loaddbNSFP {
	my $dbnsfp_annot_ref = $_[0];
	my $desiredchr = $_[1];
	my $dbnsfp_handle = $_[2];
	my $dbnsfp_keepcol_ref = $_[3];
	my $dbnsfp_nextchrline;
	
	while (<$dbnsfp_handle>) {
		$dbnsfp_nextchrline = $_;
		my @line = split("\t", $dbnsfp_nextchrline);
		my ($thischr, $pos, $ref, $alt) = @line[0..3];
		my $dbnsfp_lookup = join("_", ($thischr, $pos, $ref, $alt));
		
		if ($thischr ne $desiredchr) {
			last;
		} else {
			# parse and save desired fields
			${$dbnsfp_annot_ref}{$dbnsfp_lookup} = join("\t", @line[@{$dbnsfp_keepcol_ref}]);
		}
	}
	
	return ($dbnsfp_nextchrline, $dbnsfp_annot_ref);
}





sub loadGeneAnnot {
	my ($string, $geneannot_ref) = @_;
	my @line = split ("\t", $_);
	my ($thischr, $thisstart, $thisend) = @line[0..2];
	
	my $lookup = my $thiscoord = join('_', @line[0..5]);
	
	# 0 chromosome	begin	end	varType	reference	call	xRef
	# 7 geneId	mrnaAcc	proteinAcc	symbol	orientation	component	componentIndex
	# 14 codingRegionKnown	impact	nucleotidePos	proteinPos	annotationRefSequence
	# 19 sampleSequence	genomeRefSequence
	
	my $annot = $line[10];
	for (my $i=11; $i<=20; $i++) {
		if ($line[$i]) {
			$annot .= "\t$line[$i]";
		} else {
			$annot .= "\t";
		}
	}
	${$geneannot_ref}{$lookup} = $annot;
}
                                                                                                                                                                           





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--testvar\ttestvariants file\n";
	print "\t--genelistvar\tgenelistvar file (produced by Generate_Gene_ListVariants_From_Gene_File_0_1_4.pl)\n";
	print "\t--dbnsfp\tdbnsfp.zip file (from http://sites.google.com/site/jpopgen/dbNSFP)\n";
	print "\t--o\toutput <output.popinfo.tsv> file\n";
	die;
}