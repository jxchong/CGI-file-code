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


my ($popinfodir, $popinfofilenames, $inputfile, $outputfile);



GetOptions(
	'popinfodir=s' => \$popinfodir, 
	'popinfofilenames=s' => \$popinfofilenames, 
	'in=s' => \$inputfile,
	'out=s' => \$outputfile,
);

if (!$popinfofilenames) {
	$popinfofilenames = 'all.2012-09-20.popinfo.thischr.tsv';
}

if (!defined $popinfodir) {
	optionUsage("option --popinfodir not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $inputfile) {
	optionUsage("options --in not defined\n");
}


# my $highestchr;
my %desiredvariants;
open (INPUT, "$inputfile") or die "Cannot read $inputfile file.\n";
while ( <INPUT> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($chr, $start, $stop) = split("\t", $_);
	if ($chr =~ 'chr') {
		$chr =~ s/chr//;
	}
	$start =~ s/,//g;
	$start -= 1;
	$stop =~ s/,//g;
	push(@{$desiredvariants{$chr}}, [$start, $stop]);
}
close INPUT;



my $isFirstVariant = 1;
open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";

# iterate over each chromosome to fetch variants from each coordinate region
foreach my $chrnum (sort { $a <=> $b} keys %desiredvariants) {
	my $chrname = $chrnum;
	if ($chrname !~ 'chr') {
 		$chrname = "chr$chrnum";
 	}
	my @targetcoords = @{$desiredvariants{$chrnum}};
	# to improve efficiency, get start position of first target region and end position of last target region
	my @sortedcoords = sort { $a->[0] <=> $b->[0] } @targetcoords;
	my $firstpos = $sortedcoords[0]->[0];
	my $lastpos = $sortedcoords[$#sortedcoords]->[1];
	
	# open corresponding popinfo file for this chromosome
	my $popinfofile = $popinfofilenames;
	$popinfofile =~ s/thischr/$chrname/;
	if ("$popinfodir/$popinfofile" =~ /.bz2/) {
		open (POPINFOFILE, "bzcat $popinfodir/$popinfofile |") or die "Cannot read popinfo from $popinfodir/$popinfofile: $!\n";
	} else {
		open (POPINFOFILE, "$popinfodir/$popinfofile") or die "Cannot read popinfo from $popinfodir/$popinfofile: $!\n";
	}

	print STDERR "Fetching variants from chromosome $chrnum\n";	
	my $headerline = <POPINFOFILE>;
	$headerline =~ s/\s+$//;					# Remove line endings
	my @popinfoIDs = split("\t", $headerline);
	
	if ($isFirstVariant == 1) {
		my @header = split("\t", $headerline);
		print OUT join("\t", @header)."\n";
		$isFirstVariant = 0;
	}

	# read popinfo file line by line and check if variants are within your targetted regions
	my $matchedvariants = 0;
	while ( <POPINFOFILE> ) {
		$_ =~ s/\s+$//;
		my @line = split ("\t", $_);
		my ($thischr, $thisstart, $thisend) = @line[0..2];
		
		# skip lines that aren't in any of the targetted regions
		if ($thisstart < $firstpos) {
			next;
		} elsif ($thisend > $lastpos) {
			last;
		} else {
			# check if variant is within one of the targetted regions
			for (my $varnum=0; $varnum<=$#targetcoords; $varnum++) {
				my ($desiredstart, $desiredend, $desiredref) = @{$targetcoords[$varnum]};
				if ($thisstart >= $desiredstart && $thisend <= $desiredend) {
					$matchedvariants++;
					print OUT "$thischr\t$thisstart\t$thisend";
					for (my $i=3; $i<=$#line; $i++) {
						print OUT "\t$line[$i]";
					}
					print OUT "\n";
				}
			}
		}
	}
	close POPINFOFILE;

	print STDERR "Found $matchedvariants variants within desired coordinates on chromosome $chrnum\n";			
}

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