#!perl
#
# Description:
#
# Usage: perl untitled.pl
#
#
# Created by Jessica Chong on 2012-04-28

use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile);

GetOptions(
	'i=s' => \$inputfile, 
	'o=s' => \$outputfile,
);

open (OUT, ">$outputfile") or die "Cannot write to $outputfile.\n";
for (my $chr=1; $chr<=22; $chr++) {
	my $counter = 1;
	open (FILE, "bzcat $inputfile |") or die "Cannot read $inputfile file.\n";
	if ($chr == 1) {
		my $header = <FILE>;
		print OUT "$header";
	}
	while ( <FILE> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my @line = split ("\t", $_);
		if ($counter < 1000) {
			$counter++;
			print OUT join("\t", @line)."\n";
		} else {
			last;
		}
	}
	close FILE;
}


close OUT;




if (!defined $inputfile) {
	optionUsage("option -i not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option -o not defined\n");
} 

sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--i\tinput file\n";
	print "\t--o\toutput file\n";
	die;
}