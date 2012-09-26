#!perl
#
# Description: Add frequency, genotype class counts, and gene annotation info (if available) to an arbitrary file.
#
#
#
# Created by Jessica on 2012-04-05

use strict;
use warnings;
use Getopt::Long;

my ($inputfile, $outputfile, $popinfodir, $popinfoprefix);

GetOptions(
	'popinfodir=s' => \$popinfodir,
	'popinfoprefix=s' => \$popinfoprefix,
	'i=s' => \$inputfile, 
	'o=s' => \$outputfile,
);


if (!defined $inputfile) {
	optionUsage("option -i not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option -o not defined\n");
} elsif (!defined $popinfodir) {
	optionUsage("option -popinfodir not defined\n");
} elsif (!defined $popinfoprefix) {
	optionUsage("option -popinfoprefix not defined\n");
}


my (%allelefreqs, %genocounts);
my $allelefreqs_ref = \%allelefreqs;
my $genocounts_ref = \%genocounts;

my $currchr = 0;

open (OUT, ">$outputfile") or die "Cannot write to $outputfile.\n";
open (FILE, "$inputfile") or die "Cannot open $inputfile file.\n";
my $header = <FILE>;
$header =~ s/\s+$//;					# Remove line endings
my @headercontents = split("\t", $header);
print OUT join("\t", @headercontents)."\tMAF\tnhomref\tnhet\tnhomalt\tnhommaj\tnhommin\tCR\tgenesymbol\torientation\tcomponent\tcomponentIndex\thasCodingRegion\timpact\tnucleotidePos\tproteinPos\tannotationRefSequence\tsampleSequence\tgenomeRefSequence\n";

while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	
	print OUT join("\t", @line);
	if (@line < @headercontents) {
		my $diff = (scalar(@headercontents) - scalar(@line));
		
		for (my $i=1; $i<=$diff; $i++) {
			print OUT "\t";
		}
	}

	my $freq = 'NA';
	my @genotypes = (('NA') x 8);
	my $thiscoord = join('_', @line[0..4]);
	my $thischr = $line[0];
	
	if ($currchr ne $thischr) {
		$currchr = $thischr;
		($allelefreqs_ref, $genocounts_ref) = readPopInfo($currchr, $popinfodir, $popinfoprefix, $allelefreqs_ref, $genocounts_ref);
		print "Annotating variants on $currchr\n";
		# if ($thischr eq 'chr2') {			# DEBUG
		# 	exit;
		# }
	}
	
	if ($line[9]) {
		my $xref = $line[9];
		my @rsids = split(";", $xref);
		foreach my $rsidstring (@rsids) {
			my ($dbsnp, $rsid) = split(":", $rsidstring);
			if (exists $allelefreqs{$rsid} && $freq eq 'NA') {
				$freq = $allelefreqs{$rsid};
			}
			if (exists $genocounts{$rsid} && $genotypes[0] eq 'NA') {
				@genotypes = @{$genocounts{$rsid}};
			}
		}
	} else {
		if (exists $allelefreqs{$thiscoord}) {
			$freq = $allelefreqs{$thiscoord};
		} 
		if (exists $genocounts{$thiscoord}) {
			@genotypes = @{$genocounts{$thiscoord}};
		}
	}
	
	print OUT "\t$freq\t".join("\t", @genotypes)."\n";
}
close FILE;
close OUT;

print STDERR "done\n";








sub readPopInfo {
	my ($loadchr, $popinfodir, $popinfoprefix) = @_;
	
	%allelefreqs = ();
	%genocounts = ();
	
	print STDERR "Reading in popinfo for $loadchr\n";
	if (-e "$popinfodir/$popinfoprefix.$loadchr.bz2" ) {
		open (FREQ, "bzcat $popinfodir/$popinfoprefix.$loadchr.bz2 |") or die "Can't read annotation file $popinfodir/$popinfoprefix.$loadchr.bz2: $!\n";	
	} else {
		open (FREQ, "$popinfodir/$popinfoprefix.$loadchr.tsv") or die "Can't read annotation file $popinfodir/$popinfoprefix.$loadchr.tsv: $!\n";
	}
	my $annotationheadline = <FREQ>;
	my @annotationheader = split("\t", $annotationheadline);
	while (<FREQ>) {
		$_ =~ s/\s+$//;
		my @line = split("\t", $_);
		my $snpname = 'NA';
		if ($line[6]) {
			my $xref = $line[6];
			my @snpids = split(";", $xref);
			foreach my $snpid (@snpids) {
				my @rsids = split(":", $snpid);
				if ($rsids[0] =~ m/dbsnp/i) {
					$snpname = $rsids[1];
				}
			}
			if ($snpname eq 'NA') {
				$snpname = $line[7];
			}
		} else {
			$snpname = $line[7];
		}

		$allelefreqs{$snpname} = $line[8];

		my $geneannot = (("\t") x 11);
		if ($line[16]) {
			$geneannot = $line[16];
			for (my $i=17; $i<=20; $i++) {
				if ($line[$i]) {
					$geneannot .= "\t$line[$i]";
				} else {
					$geneannot .= "\t";
				}
			}
		}

		@{$genocounts{$snpname}} = (@line[9..11,13..15], $geneannot);
		
	}
	close FREQ;
	
	return (\%allelefreqs, \%genocounts);
}










sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--i\tinput file (produced by getvar_mastervarfilesbz2_bygeno.pl)\n";
	print "\t--popinfo\t.popinfo.tsv or .popinfo.tsv.bz2 file produced by calc_popinfo.pl\n";
	print "\t--o\toutput file\n";
	die;
}