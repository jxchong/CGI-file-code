#!perl
#
# Description:e
#
# Usage: perl untitled.pl
#
#
# Created by Jessica Chong on 2012-04-28

use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile, $column);

GetOptions(
	'i=s' => \$inputfile, 
	'o=s' => \$outputfile,
	'column=s' => \$column,
);

my $counter = 0;
my $oldchr = 'chr1';
open (OUT, ">$outputfile") or die "Cannot write to $outputfile.\n";
open (FILE, "bzcat $inputfile |") or die "Cannot read $inputfile file.\n";
my $header = <FILE>;
print OUT "$header";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	if ($line[$column] ne $oldchr) {
		$counter = 1;
		$oldchr = $line[$column];
		print OUT join("\t", @line)."\n";
	} elsif ($line[$column] eq $oldchr && $counter < 100) {
		$counter++;
		print OUT join("\t", @line)."\n";
	}
}
close FILE;
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