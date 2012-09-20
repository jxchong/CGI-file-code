#!/usr/bin/env perl
#
# Description: Make submission script for genelistvariants
#
#
#
# Created by Jessica Chong on 2012-09-20.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile);

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} 



open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "#!/bin/bash\n";
print OUT '#$-cwd'."\n\n";
print OUT "perl ~/bin/Generate_Gene_ListVariants_From_Gene_File_0_1_4.pl ";
open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my $file = $line[0];
	$file =~ s/masterVarBeta/gene/;
	print OUT "--Input_File $file "
}
close FILE;
print OUT "--Output_Dir /clustb/jxchong/CGI_WGS --Output_File all.2012-09-20.genelistvar.tsv\n\n";
print OUT "cp /clustb/jxchong/CGI_WGS/all.2012-09-20.genelistvar.tsv.bz2 /clusta/jxchong/CGI_WGS/all.2012-09-20.genelistvar.tsv.bz2\n";
close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file (findiv_mastervarpath)\n";
	print "\t--out\toutput file (output for submission script)\n";
	die;
}