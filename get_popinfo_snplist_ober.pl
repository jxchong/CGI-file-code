#!/usr/bin/env perl
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
	# 'popinfodir=s' => \$popinfodir, 
	# 'popinfofilenames=s' => \$popinfofilenames, 
	'in=s' => \$inputfile,
	'out=s' => \$outputfile,
);


if (-d "/ober_resources/CGI_98Hutterites_WGS/popinfo-2012-11-12") {
	$popinfodir = '/clusta/jxchong/CGI_WGS/popinfo-2012-11-12';
} elsif (-d "/clusta/jxchong/CGI_WGS/popinfo-2012-11-12") {
	$popinfodir = '/clusta/jxchong/CGI_WGS/popinfo-2012-11-12';
} else {
	print STDERR "Cannot access popinfo files in either directory: /ober_resources/CGI_98Hutterites_WGS/popinfo-2012-11-12 or /clusta/jxchong/CGI_WGS/popinfo-2012-11-12\n";
	die;
}

$popinfofilenames = 'all.2012-11-12.popinfo.tsv.gz';


if (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $inputfile) {
	optionUsage("options --in not defined\n");
}


print STDOUT "Note: These results (genotype counts and frequencies) use data from only 96 Hutterite genomes (not from 98) because including 2 affected children biases the MAF calculations\n\n";


# get header
open (POPINFOFILE, "zcat $popinfodir/$popinfofilenames |") or die "Cannot read popinfo from $popinfodir/$popinfofilenames: $!\n";
my $header = <POPINFOFILE>;
$header =~ s/\s+$//;					# Remove line endings
close POPINFOFILE;



open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "$header\n";
open (INPUT, "$inputfile") or die "Cannot read $inputfile file.\n";
while ( <INPUT> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($chr, $start, $stop) = split("\t", $_);
	if ($chr !~ 'chr') {
		$chr = "chr$chr";
	}
	$start =~ s/,//g;
	$start -= 1;
	$stop =~ s/,//g;
	
	print STDOUT "Fetching variants from chromosome $chr:$start-$stop\n";
	
	my $matchedvariants = 0;
	my @variants = `tabix $popinfodir/$popinfofilenames $chr:$start-$stop`;
	foreach (@variants) {
		print OUT "$_";
		$matchedvariants++;
	}
	print STDOUT "Found $matchedvariants variants between coordinates $chr:$start-$stop\n";			
	
}
close INPUT;

close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--out\toutput file\n";
	print "\t--in\tinput file, tab-delimited, coordinates assumed to be 1-based: chr, start, stop\n";
	die;
}