#!perl
#
# Description:
#
# Usage: perl untitled.pl
#
#
# Created by Jessica Chong on 2012-06-11

use strict;
use warnings;
use Getopt::Long;


my ($testvarfile, $desiredchr, $desiredstart, $desiredend, $outputfile);

GetOptions(
	'testvar=s' => \$testvarfile, 
	'chr=s' => \$desiredchr,
	'start=i' => \$desiredstart,
	'end=i' => \$desiredend,
	'out=s' => \$outputfile,
);

if (!defined $testvarfile) {
	optionUsage("option --testvar not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $desiredchr) {
	optionUsage("options --chr not defined\n");
} elsif (!defined $desiredstart) {
	optionUsage("options --start not defined\n");
} elsif (!defined $desiredend) {
	optionUsage("options --end not defined\n");
}

if ($desiredchr !~ 'chr') {
	$desiredchr = "chr$desiredchr";
}
if ($desiredend < $desiredstart) {
	print STDERR "End $desiredend is before start $desiredstart\n";
	die;
} elsif ($desiredstart > $desiredend) {
	print STDERR "Start $desiredstart is after end $desiredend\n";
	die;	
}



my (%findiv2CGI, %CGI2findiv);
open (FILE, "/ober_resources/CGI_98Hutterites_WGS/README.assembly_sample_subject.csv") or die "Cannot read README.assembly_sample_subject.csv file.\n";
<FILE>;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split (",", $_);
	$findiv2CGI{$line[2]} = $line[0];
	$CGI2findiv{$line[0]} = $line[2];
}
close FILE;


# my $desiredCGIid = $findiv2CGI{$desiredfindiv};

open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "variantId\tChr\tStart\tStop\tFindiv\tGenotype\n";
open (FILE, "bzcat $testvarfile |") or die "Cannot read $testvarfile file: $!.\n";
my $headerline = <FILE>;
$headerline =~ s/\s+$//;					# Remove line endings
my @testvarIDs = split("\t", $headerline);

my $currchr = 'NA';
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($thisvarnum, $thischr, $thisstart, $thisend) = @line[0..3];
	
	if ($currchr ne $thischr) {
		print STDERR "Reading chromosome $thischr\n";
		$currchr = $thischr;
	}
	
	if ($thischr ne $desiredchr) {
		next;
	} elsif ($thischr eq $desiredchr && $thisstart >= $desiredstart && $thisend <= $desiredend) {
		for (my $i=8; $i<=$#line; $i++) {
			print OUT "$thisvarnum\t$thischr\t$thisstart\t$thisend\t$CGI2findiv{$testvarIDs[$i]}\t$line[$i]\n";
		}
	} elsif ($thischr eq $desiredchr && $thisend > $desiredend) {
		last;
	}
}
close FILE;
close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--testvar\ttestvar file\n";
	print "\t--out\toutput file\n";
	print "\t--chr\tdesired chr, must be 'chrX, chr1' style\n";
	print "\t--start\tdesired start\n";
	print "\t--end\tdesired end\n";
	die;
}