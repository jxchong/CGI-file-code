#!perl
#
# Description: Calculate allele frequency and genotype counts from testvariants data file.
#
#
#
# Created by Jessica on 2012-04-09

use strict;
use warnings;
use Getopt::Long;

my ($inputfile, $outputfile);

GetOptions(
	'i=s' => \$inputfile, 
	'o=s' => \$outputfile,
);

if (!defined $inputfile) {
	optionUsage("option -i not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option -o not defined\n");
} 


my $currchr = 'chr1';

open (OUT, ">$outputfile.popinfo.tsv");
print OUT "chr\tbegin\tend\tvartype\tref\talt\txref\tlocusname\taltfreq\tMAF\tnhomref\tnhet\tnhomalt\tnmiss\tnhommaj\tnhommin\tCR\n";

open (FILE, "$inputfile") or die "Cannot read $inputfile file.\n";
<FILE>;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($thischr, $thisstart, $thisend) = @line[1..3];
	if ($thischr ne $currchr) {
		$currchr = $thischr;
		print "Processing chr $currchr\n";
	}
	my $thiscoord = $thischr.'_'.$thisstart.'_'.$thisend;
	my @genotypes = @line[8..$#line];
	my @genocounts = (0,0,0,0);
	foreach my $geno (@genotypes) {
		if ($thischr =~ 'chrX') {
			last;
			if ($geno eq '00')  { $genocounts[0]++; }
			if ($geno eq '01' || $geno eq '10')  { $genocounts[1]++; }
			if ($geno eq '11')  { $genocounts[2]++; }			
			if ($geno =~ 'N')  { $genocounts[3]++; }
		} elsif ($thischr =~ 'chrM' || $thischr =~ 'chrY') {
			last;
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
	
		my ($maf, $genocounts_relhutt_maj, $genocounts_relhutt_min);		# count number of each genotype class using the hutterite minor allele as the alternative allele
		if ($genocounts[0] >= $genocounts[2]) {			# more ref alleles than alt alleles
			$maf = ($genocounts[1]+2*$genocounts[2])/(2*$ngenosubj);
			$genocounts_relhutt_maj = $genocounts[0];
			$genocounts_relhutt_min = $genocounts[2];
		} else {																	# same or more alt alleles as ref alleles
			$maf = ($genocounts[1]+2*$genocounts[0])/(2*$ngenosubj);
			$genocounts_relhutt_maj = $genocounts[2];
			$genocounts_relhutt_min = $genocounts[0];
		}
	
		
		for (my $i=1; $i<=7; $i++) {
			if ($line[$i]) {
				print OUT "$line[$i]\t";
			} else {
				print OUT "\t";
			}
		}
		print OUT "$thiscoord\t$altfreq\t$maf\t".join("\t", @genocounts)."\t$genocounts_relhutt_maj\t$genocounts_relhutt_min\t$callrate\n";		
	}

	
}
close FILE;
close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--i\tinput testvariants file\n";
	print "\t--o\toutput <output>.popinfo.tsv file\n";
	die;
}