#!perl
#
# Description: Get all fully called, non-complex variants from a mastervar file matching a particular genotype
# 	Can get variants within a given genomic coordinate or throughout the entire genome
#
#
#
# Created by Jessica on 2012-03-26

use strict;
use warnings;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);
use Getopt::Long;


my ($mastervarfile, $desiredgeno, $targetchr, $targetstart, $targetend, $outputfile);

GetOptions(
	'i=s' => \$mastervarfile, 
	'genotype=s' => \$desiredgeno,
	'chr:s' => \$targetchr,
	'start:i' => \$targetstart,
	'end:i' => \$targetend,
	'o=s' => \$outputfile,
);


if (!defined $mastervarfile) {
	optionUsage("option -i not defined\n");
} elsif (!defined $desiredgeno) {
	optionUsage("option -genotype not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option -o not defined\n");
} elsif (!defined $targetchr && (!defined $targetstart || !defined $targetend)) {
	optionUsage("--chr provided but start or stop not defined\n");
}

sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--i\tinput masterVarBeta-ASM* file\n";
	print "\t--genotype\tgenotype desired, either: homref / het / homalt / all\n";
	print "\t--chr\tchromosome (optional)\n";
	print "\t--start\tstart position (optional)\n";
	print "\t--end\tend position (optional)\n";
	print "\t--o\toutput file\n";
	die;
}

if ($mastervarfile !~ m/masterVarBeta/) {
	print STDERR "input file ($mastervarfile) is not a mastervar file\n";
	die;
}
if ($targetchr !~ m/chr/) {
	$targetchr = 'chr'.$targetchr;
}


my $targetchrnum = $targetchr;
$targetchrnum =~ s/chr//;

my %allowedvartypes = (
	'snp' => 1,
	'ins' => 1,
	'del' => 1,
	'sub' => 1,
	'ref' => 1,
);

my %allowedzygosity = (
	'hap' => 1,
	'hom' => 1,
	'het-ref' => 1,
	'het-alt' => 1,
);


my $currchr = 'NA';
my @keepfields = (2..9, 18, 19, 25..28);
open (OUT, ">$outputfile");
# my $bzobj = new IO::Uncompress::Bunzip2 $mastervarfile or die "IO::Uncompress::Bunzip2 failed: $Bunzip2Error\n";
open (BZ, "bzcat $mastervarfile |") or die "Cannot read mastervarfile $mastervarfile: $!\n";
# while (<$bzobj>) {
while (<BZ>) {	
	my $nextline = $_;
	if ($nextline =~ m/^#/ || $nextline =~ m/^\s*$/) {					# skip header lines
		next;
	} elsif ($_ =~ m/^\>/) {				
		$nextline =~ s/\s+$//;									# Remove line endings
		my @headerline = split("\t", $nextline);
		print OUT join("\t", @headerline[@keepfields])."\n";
	} else {
		$nextline =~ s/\s+$//;										# Remove line endings
		my @line = split ("\t", $nextline);
		my ($thischr, $thisstart, $thisend, $zygosity, $vartype, $refallele, $a1, $a2) = @line[2..9];
		
		if ($currchr ne $thischr) {
			print STDERR "At $thischr: looking for $targetchr:$targetstart-$targetend\n";
			$currchr = $thischr;
		}
		
		if (defined $targetchr && defined $targetstart && defined $targetend) {
			my $thischrnum = $thischr;
			$thischrnum =~ s/chr//;

			if ($thischrnum > $targetchrnum || (($thischr eq $targetchr ) && ($thisend > $targetend))) {
				last;
			}

			if ($thischr eq $targetchr && $thisstart>=$targetstart && $thisend<=$targetend) {
				if (exists $allowedvartypes{$vartype} && exists $allowedzygosity{$zygosity} && testGenoClass($desiredgeno, $zygosity, $vartype) == 1) {
					for (my $i=0; $i<=$#keepfields; $i++) {
						my $fieldnum = $keepfields[$i];
						if (defined $line[$fieldnum]) {
							print OUT "$line[$fieldnum]\t";
						} else {
							print OUT "\t";
						}
					}
					print OUT "\n";
				}
			}
		} else {			# no chr coordinates specified, so look at all variants in genome
			if (exists $allowedvartypes{$vartype} && exists $allowedzygosity{$zygosity} && testGenoClass($desiredgeno, $zygosity, $vartype) == 1) {
				for (my $i=0; $i<=$#keepfields; $i++) {
					my $fieldnum = $keepfields[$i];
					if (defined $line[$fieldnum]) {
						print OUT "$line[$fieldnum]\t";			# rs11934749
					} else {
						print OUT "\t";
					}
				}
				print OUT "\n";
			}
		}
	}
}

# $bzobj->close();
close BZ;

close OUT;



print STDERR "done\n";


sub testGenoClass {
	my ($testtype, $zygosity, $vartype) = @_;
	my $testresult;
	
	if ($testtype eq 'het') {
		if ($zygosity eq 'het-ref' && $vartype ne 'ref') {
			$testresult = 1;
		} else {
			$testresult = 0;
		}
	}
	
	if ($testtype eq 'homalt') {
		if (($zygosity eq 'hom' || $zygosity eq 'het-alt') && $vartype ne 'ref') {
			$testresult = 1;
		} else {
			$testresult = 0;
		}
	}
	
	if ($testtype eq 'homref') {
		if ($zygosity eq 'hom' && $vartype eq 'ref') {
			$testresult = 1;
		} else {
			$testresult = 0;
		}
	}
	
	if ($testtype eq 'all') {
		$testresult = 1;
	}
	
	return $testresult;
}

