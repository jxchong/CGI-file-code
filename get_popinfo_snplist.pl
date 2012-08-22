#!perl
#
# Description: Get all variants within certain coordinates from a popinfo file.
#
#
#
# Created by Jessica X Chong on 2012-06-19


use strict;
use warnings;
use Getopt::Long;


my ($popinfofile, $inputfile, $outputfile);

GetOptions(
	'popinfo=s' => \$popinfofile, 
	'in=s' => \$inputfile,
	'out=s' => \$outputfile,
);



if (!defined $popinfofile) {
	optionUsage("option --popinfo not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $inputfile) {
	optionUsage("options --in not defined\n");
}


# my $highestchr;
my %desiredvariants;
open (FILE, "$inputfile") or die "Cannot read $inputfile file.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($chr, $start, $stop) = split("\t", $_);
	if ($chr !~ 'chr') {
		$chr = "chr$chr";
	}
	$start -= 1;
	push(@{$desiredvariants{$chr}}, [$start, $stop]);

	# my $thischrnum = $chr;
	# $thischrnum =~ s/chr//;
	# if ($thischrnum > $highestchr) {
	# 	$highestchr = $thischrnum;
	# }
}
close FILE;






open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (FILE, "bzcat $popinfofile |") or die "Cannot read $popinfofile file: $!.\n";
my $headerline = <FILE>;
$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("\t", $headerline);
my @popinfoIDs = split("\t", $headerline);
my @thischrvariants;

print OUT join("\t", @header)."\n";
my $currchr = 'NA';
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($thischr, $thisstart, $thisend) = @line[0..2];
	
	if ($currchr ne $thischr) {
		print STDERR "Reading chromosome $thischr\n";
		$currchr = $thischr;
		if (exists $desiredvariants{$thischr}) {
			@thischrvariants = @{$desiredvariants{$thischr}};
			print STDERR "There are ".scalar(@thischrvariants)." desired variants on $currchr\n";			
		}
	}
	
	if (!exists $desiredvariants{$thischr}) {
		next;
	} else {
		for (my $varnum=0; $varnum<=$#thischrvariants; $varnum++) {
			my ($desiredstart, $desiredend) = @{$thischrvariants[$varnum]};
			if ($thisstart >= $desiredstart && $thisend <= $desiredend) {
				print OUT "$thischr\t$thisstart\t$thisend";
				for (my $i=3; $i<=$#line; $i++) {
					print OUT "\t$line[$i]";
				}
				print OUT "\n";
			}
		}
		
	} 
	
	
}
close FILE;
close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--popinfo\tpopinfo file\n";
	print "\t--out\toutput file\n";
	print "\t--in\tinput file, tab-delimited, coordinates assumed to be 1-based: chr, start, stop\n";
	die;
}