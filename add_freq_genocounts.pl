#!perl
#
# Description: Add frequency and genotype class counts.
#
#
#
# Created by Jessica on 2012-04-05

use strict;
use warnings;
use Getopt::Long;

my ($inputfile, $outputfile, $popinfofile);

GetOptions(
	'popinfo=s' => \$popinfofile,
	'i=s' => \$inputfile, 
	'o=s' => \$outputfile,
);


if (!defined $inputfile) {
	optionUsage("option -i not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option -o not defined\n");
} elsif (!defined $popinfofile) {
	optionUsage("option -popinfo not defined\n");
} 


print STDERR "Reading in allele freqs\n";
my (%allelefreqs, %genocounts);
open (FREQ, "$popinfofile") or die "Can't read freqfile $popinfofile\n";
<FREQ>;
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
	
	$allelefreqs{$snpname} = $line[9];
	@{$genocounts{$snpname}} = @line[10..12,14..16];
	
	# if ($snpname eq 'chr1_95621744_95621747') {
	# 	# for (my $i=0; $i<=$#line; $i++) {
	# 	# 	print "$i = $line[$i]\n";
	# 	# }
	# 	# print "freq: $allelefreqs{$snpname}\n";
	# 	# print "genocounts: @{$genocounts{$snpname}}\n";
	# 	last;
	# }
}
close FREQ;





open (OUT, ">$outputfile") or die "Cannot write to $outputfile.\n";
open (FILE, "$inputfile") or die "Cannot open $inputfile file.\n";
my $header = <FILE>;
$header =~ s/\s+$//;					# Remove line endings
my @headercontents = split("\t", $header);
print OUT join("\t", @headercontents)."\tMAF\tnhomref\tnhet\tnhomalt\tnhommaj\tnhommin\tCR\n";

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
	my @genotypes = (('NA') x 6);
	my $thiscoord = $line[0].'_'.$line[1].'_'.$line[2];
	
	if ($line[8]) {
		my $xref = $line[8];
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
	
	# if ($thiscoord eq 'chr1_95621744_95621747') {
	# 	print "freq: $allelefreqs{$thiscoord}\n";
	# 	@genotypes = @{$genocounts{$thiscoord}};
	# 	print "genocounts: @{$genocounts{$thiscoord}}\n";
	# 	print "genotypes: @genotypes\n";
	# 	exit;
	# }
}
close FILE;
close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--i\tinput file (produced by getvar_mastervarfilesbz2_bygeno.pl)\n";
	print "\t--popinfo\t.popinfo.tsv file produced by calc_popinfo.pl\n";
	print "\t--o\toutput file\n";
	die;
}