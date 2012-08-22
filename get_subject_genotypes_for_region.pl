#!perl
#
# Description: From a given testvar file, get all variants in a specific region and produce a PLINK-like output of genotypes
#
#
#
# Created by Jessica X Chong on 2012-06-11

use strict;
use warnings;
use Getopt::Long;


my ($testvarfile, $desiredchr, $desiredstart, $desiredend, $outputfile, $translateIDfile);

GetOptions(
	'testvar=s' => \$testvarfile, 
	'chr=s' => \$desiredchr,
	'start=i' => \$desiredstart,
	'end=i' => \$desiredend,
	'out=s' => \$outputfile,
	'idtranslate=s' => \$translateIDfile,
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
} elsif (!defined $translateIDfile) {
	optionUsage("options --idtranslate not defined")
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



my (%subjectID2CGI, %CGI2subjectID);
my $translateIDfile; 
open (FILE, "$translateIDfile") or die "Cannot read $translateIDfile file.\n";
<FILE>;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split (",", $_);
	$subjectID2CGI{$line[2]} = $line[0];
	$CGI2subjectID{$line[0]} = $line[2];
}
close FILE;


# my $desiredCGIid = $subjectID2CGI{$desiredsubjectID};

open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "Chr\tStart\tStop\tsubjectID\tGenotype\n";
open (FILE, "bzcat $testvarfile |") or die "Cannot read $testvarfile file: $!.\n";
my $headerline = <FILE>;
$headerline =~ s/\s+$//;					# Remove line endings
my @testvarIDs = split("\t", $headerline);

my $currchr = 'NA';
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($thischr, $thisstart, $thisend) = @line[1..3];
	
	if ($currchr ne $thischr) {
		print STDERR "Reading chromosome $thischr\n";
		$currchr = $thischr;
	}
	
	if ($thischr ne $desiredchr) {
		next;
	} elsif ($thischr eq $desiredchr && $thisstart >= $desiredstart && $thisend <= $desiredend) {
		for (my $i=8; $i<=$#line; $i++) {
			print OUT "$thischr\t$thisstart\t$thisend\t$CGI2subjectID{$testvarIDs[$i]}\t$line[$i]\n";
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
	print "\t--idtranslate\t file containing translation between CGI and your IDs\n";
	die;
}