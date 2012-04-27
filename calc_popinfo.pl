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


my ($inputfile, $outputfile, $geneannotationlistvar, %geneannotations);

GetOptions(
	'testvar=s' => \$inputfile, 
	'o=s' => \$outputfile,
	'genelistvar=s' => \$geneannotationlistvar,
);

if (!defined $inputfile) {
	optionUsage("option --testvar not defined\n");
} elsif (!defined $geneannotationlistvar) {
	optionUsage("option --genelistvar not defined\n");
}	elsif (!defined $outputfile) {
	optionUsage("option --o not defined\n");
}

# load in gene annotations from files
print STDERR "Preloading gene annotations from $geneannotationlistvar\n";
open (GENELISTVAR, "$geneannotationlistvar") or die "Cannot read $geneannotationlistvar\n";
<GENELISTVAR>;
while (<GENELISTVAR>) {
	$_ =~ s/\s+$//;
	loadGeneAnnot($_, \%geneannotations);
}
close GENELISTVAR;


open (OUT, ">$outputfile");
print OUT "chr\tbegin\tend\tvartype\tref\talt\txref\taltfreq\tMAF\tnhomref\tnhet\tnhomalt\tnmiss\tnhommaj\tnhommin\tCR\tgenesymbol\torientation\tcomponent\tcomponentIndex\thasCodingRegion\timpact\tnucleotidePos\tproteinPos\tannotationRefSequence\tsampleSequence\tgenomeRefSequence\n";
open (FILE, "$inputfile") or die "Cannot read $inputfile file.\n";
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
	
	if ($thischr ne $currchr) {
		print STDERR "Working on $thischr\n";			
		$currchr = $thischr;
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
		print OUT "$altfreq\t$maf\t".join("\t", @genocounts)."\t$genocounts_relhutt_maj\t$genocounts_relhutt_min\t$callrate\t$geneannot\n";		
	}
}
close FILE;
close OUT;





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
	print "\t--o\toutput <output.popinfo.tsv> file\n";
	die;
}