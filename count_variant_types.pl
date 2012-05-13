#!perl
#
# Description: Count number of variants in popinfo file (made by calc_popinfo.pl) that are novel, in gene region, nonsense, missense, etc.
#
#
#
# Created by Jessica Chong on 2012-05-11

use strict;
use warnings;
use Getopt::Long;


my ($popinfofile, $outputfile);

GetOptions(
	'popinfo=s' => \$popinfofile, 
	'o=s' => \$outputfile,
);

my $xrefcol = 6;
my $callratecol = 15;
my $componentcol = 18;
my $impactcol = 21;

my $highCR = 0;
my $novelcount = 0;
my %component;
my %impact;

open (FILE, "bzcat $popinfofile |") or die "Cannot read $popinfofile file.\n";
my $headerline = <FILE>;
my @header = split("\t", $headerline);
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	if ($line[$callratecol] >= 0.9) {
		$highCR++;
		
		if (!$line[$xrefcol]) {
			$novelcount++;
		}
		
		if ($line[$componentcol]) {
			$component{$line[$componentcol]} += 1;
		}
		
		if ($line[$impactcol]) {
			$impact{$line[$impactcol]} += 1;
		}
		
	}
}
close FILE;

open (OUT, ">$outputfile") or die "Cannot write to $outputfile.\n";
print OUT "$highCR with CR>=90%\n";
print OUT "$novelcount with CR>=90% and novel\n";
print OUT "COMPONENT:\n";
while (my($type, $count) = each %component) {
	print OUT "$type\t$count\n";
}
print OUT "IMPACT:\n";
while (my($type, $count) = each %impact) {
	print OUT "$type\t$count\n";
}
close OUT;




if (!defined $popinfofile) {
	optionUsage("option -popinfo not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option -o not defined\n");
} 

sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--popinfo\tpopinfo file\n";
	print "\t--o\toutput file\n";
	die;
}