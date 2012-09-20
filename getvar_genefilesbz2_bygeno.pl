#!perl
#
# Description:
#
# Usage: perl untitled.pl
#
#
# Created by Jessica on 2012-03-26

use strict;
use warnings;
# use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);
use Getopt::Long;


my ($genevarfile, $desiredgeno, $targetchr, $targetstart, $targetend, $outputfile);

GetOptions(
	'i=s' => \$genevarfile, 
	'genotype=s' => \$desiredgeno,
	'chr:s' => \$targetchr,
	'start:i' => \$targetstart,
	'end:i' => \$targetend,
	'o=s' => \$outputfile,
);


if (!defined $genevarfile) {
	optionUsage("option --i not defined\n");
} elsif (!defined $desiredgeno) {
	optionUsage("option --genotype not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --o not defined\n");
} 

sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--i\tinput gene-ASM* file\n";
	print "\t--genotype\tgenotype desired, either: homref / het / homalt\n";
	print "\t--chr\tchromosome (optional)\n";
	print "\t--start\tstart position (optional)\n";
	print "\t--end\tend position (optional)\n";
	print "\t--o\toutput file\n";
	die;
}


my $header;
my $prevallele;

my $countlines = 0;
open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (IN, "bzcat $genevarfile |") or die "Cannot read $genevarfile: $!.\n";
while (<IN>) {
	$countlines++;
	my $nextline = $_;
	if ($nextline =~ m/^#/ || $nextline =~ m/^\s*$/) {					# skip header lines
		next;
	} elsif ($_ =~ m/^\>/) {				
		$nextline =~ s/\s+$//;									# Remove line endings
		my @headerline = split("\t", $nextline);
		print OUT join("\t", @headerline[3..5]);
		print OUT "\tref\ta1\ta2\t";
		print OUT join("\t", @headerline[9,13..23])."\taminochange\n";
	} else {
		$nextline =~ s/\s+$//;										# Remove line endings
		my @line = split ("\t", $nextline);
		my ($a1, $a2);
		my ($allelenum, $thischr, $thisstart, $thisend, $vartype, $refallele, $calledallele) = @line[2..8];
		my $varloc = $line[15];
		
		# if ($countlines > 50) {			 # DEBUG
		# 	print STDERR $nextline;
		# 	print STDERR "$allelenum\n$thischr\n$thisstart\n$thisend\n$vartype\n$refallele\n$calledallele\n";
		# 	exit;
		# }
		
		if (!defined $targetchr || ($thischr eq $targetchr && $thisstart>=$targetstart && $thisend<=$targetend)) {
			if ($thischr eq 'chrX' || $thischr eq 'chrY' || $thischr eq 'chrM') {
				last;
			} elsif ($vartype eq 'snp' || $vartype eq 'ref') {
				if ($allelenum == 1) {
					$prevallele = $calledallele;
					# print STDERR "for $allelenum, call=$calledallele\n";		 # DEBUG
				} elsif ($allelenum == 2 && testGenoClass($desiredgeno, $prevallele, $calledallele, $refallele) == 1) {
					print OUT "$thischr\t$thisstart\t$thisend\t$refallele\t$prevallele\t$calledallele\t";
					if (defined $line[9]) {
						print OUT "$line[9]\t";
					} else {
						print OUT "\t";
					}
					for (my $i=13; $i<=23; $i++) {
						if (defined $line[$i]) {
							print OUT "$line[$i]\t";
						} else {
							print OUT "\t";
						}
					}
					if (defined $line[20] && defined $line[22] && defined $line[23]) {
						print OUT "$line[23]$line[20]$line[22]";
					} else {
						print OUT "\t";
					}
					print OUT "\n";
				}		
			}
		}
	}
}

close IN;
close OUT;


print STDERR "done\n";


sub testGenoClass {
	my ($testtype, $prevallele, $calledallele, $refallele) = @_;
	my $testresult;
	
	if ("$prevallele$calledallele" =~ '0' || "$prevallele$calledallele" =~ 'N') {
		$testresult = '0';
	}	
	if ($testtype eq 'het') {
		$testresult = isHet($prevallele, $calledallele, $refallele);
	}
	if ($testtype eq 'homalt') {
		$testresult = isHomAlt($prevallele, $calledallele, $refallele);
	}
	if ($testtype eq 'homref') {
		$testresult = isHomRef($prevallele, $calledallele, $refallele);
	}
	
	return $testresult;
}


sub isHet {
	my ($a1, $a2, $ref) = @_;
	
	if ($a1 ne $a2 && "$a1$a2" =~ $ref) {
		return 1;
	} else {
		return 0;
	}
}

sub isHomAlt {
	my ($a1, $a2, $ref) = @_;
	
	if ($a1 eq $a2 && $a1 ne $ref) {
		return 1;
	} else {
		return 0;
	}
}

sub isHomRef {
	my ($a1, $a2, $ref) = @_;
	
	if ($a1 eq $a2 && $a1 eq $ref) {
		return 1;
	} else {
		return 0;
	}
}

