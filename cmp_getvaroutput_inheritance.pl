#!perl
#
# Description: Compare variants by chr and start position in files outputted by getvar_mastervarfilesbz2_bygeno.pl
#
#
#
# Created by Jessica on 2012-03-27

use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile);

GetOptions(
	'i=s' => \$inputfile, 
	'o=s' => \$outputfile,
);


if (!defined $inputfile) {
	optionUsage("option -i not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option -o not defined\n");
}

sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--i\tinput file: list of files to compare; first three columns must be chr, start, stop\n";
	print "\t--o\toutput file\n";
	die;
}



my @files;
open (IN, "$inputfile");
while (<IN>) {
	$_ =~ s/\s+$//;
	if ($_ !~ '#') {
		push(@files, $_);
	}
}
close IN;


my %varpresent;
for (my $i=0; $i<=$#files; $i++) {
	print STDERR "Looking at variants in $files[$i]\n";
	my $countvar = 0;
	open (FILE, "$files[$i]") or die "Cannot open $files[$i].\n";
	<FILE>;
	while ( <FILE> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my @line = split ("\t", $_);
		# my $pos = $line[0].'_'.$line[1].'_'.$line[2];
		my $lookup = join('_', @line[0..2,5]);
		${$varpresent{$lookup}}[$i] += 1;
		for (my $j=0; $j<=$#files; $j++) {
			${$varpresent{$lookup}}[$j] += 0;
		}
		$countvar++;
		if ($line[0] eq 'chr5') {
			last;
		}
	}
	close FILE;	
	print STDERR ".. checked $countvar variants\n";
}

# while (my($pos,$countref) = each %varpresent) {
# 	print STDERR "variant at $pos, @{$countref}\n";
# }

open (OUT, ">$outputfile") or die "Cannot write to $outputfile\n";
my $countvarmatches = 0;
open (FILE, "$files[$#files]") or die "Cannot open $files[$#files].\n";
my $header = <FILE>;
print OUT $header;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	# my $pos = $line[0].'_'.$line[1].'_'.$line[2];
	my $lookup = join('_', @line[0..2,5]);
	if ($line[0] eq 'chr5') {
		last;
	}
	if (exists $varpresent{$lookup}) {
		my $countmatches = 0;
		for (my $i=0; $i<=$#files; $i++) {
			if (${$varpresent{$lookup}}[$i] >= 1) {
				$countmatches++;
			}
		}
		
		if ($countmatches == scalar(@files)){
			print OUT join("\t", @line)."\n";
			$countvarmatches++;
		}
	}
}
close FILE;
close OUT;
print STDERR "Found $countvarmatches variants in common\n";
