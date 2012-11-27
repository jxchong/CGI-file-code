#!/usr/bin/env perl
#
# Description: Split output by chromosome (first column contents)
#
#
#
# Created by Jessica Chong on 2012-09-20.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputprefix, $chrcol);

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputprefix,
	'chrcol=i' => \$chrcol,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputprefix) {
	optionUsage("option --out not defined; format=xxxxx.chr.tsv where chr will be replaced with chromosome number\n");
} elsif (!defined $chrcol) {
	optionUsage("option --chrcol not defined; indicate column number (1-based) with the chromosome number\n");
}


my @filenames;

my $currchr = 0;

if ($inputfile =~ /.bz2$/) {
	open (FILE, "bzcat $inputfile |") or die "Cannot read $inputfile: $!\n";
} elsif ($inputfile =~ /.gz$/) {
	open (FILE, "zcat $inputfile |") or die "Cannot read $inputfile: $!\n";
} else {
	open (FILE, "$inputfile") or die "Cannot read $inputfile: $!\n";
}
my $header = <FILE>;
$header =~ s/\s+$//;					# Remove line endings
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my $thischr = $line[($chrcol-1)];
	if ($currchr ne $thischr) {
		if ($currchr ne '0') {
			close OUT;	
		}
		$currchr = $thischr;
		
		my $outputfile = $outputprefix;
		$outputfile =~ s/chr/$thischr/;
		push(@filenames, $outputfile);
		print STDERR "Working on chromosome $thischr\n";
		open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
		print OUT "$header\n";
	}
	
	print OUT join("\t", @line)."\n";
}
close FILE;


print STDERR "Compressing files\n";

foreach my $filename (@filenames) {
	`~/bin/bgzip $filename`;
}





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	die;
}