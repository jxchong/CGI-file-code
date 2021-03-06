#!perl
#
# Description: Calculate allele frequency and genotype counts from testvariants data file.  Combine with gene annotation info (if available) into one file.  Everything is much faster if it's compressed in .bz2 format!!
#
# Uses dbNSFP http://sites.google.com/site/jpopgen/dbNSFP
#
# Created by Jessica X Chong on 2012-04-09

use strict;
use warnings;
use Getopt::Long;


my ($testvarfile, $outputfile, $geneannotationlistvar, %geneannotations, $dbnsfpfile);
my %dbnsfp_annots;
my $dbnsfp_annot_ref = \%dbnsfp_annots;
my @dbnsfp_keepcol = (14, 17, 18, 21..34, 37, 47, 48);		# choose columns to keep from dbNSFP

my $countSNPnogenosubj = 0;

GetOptions(
	'testvar=s' => \$testvarfile, 
	'o=s' => \$outputfile,
	'genelistvar=s' => \$geneannotationlistvar,
	'dbnsfp=s' => \$dbnsfpfile,
);

if (!defined $testvarfile) {
	optionUsage("option --testvar not defined\n");
} elsif (!defined $geneannotationlistvar) {
	optionUsage("option --genelistvar not defined\n");
}	elsif (!defined $outputfile) {
	optionUsage("option --o not defined\n");
} elsif (!defined $dbnsfpfile) {
	optionUsage("option --dbnsfp file not defined\n");
} 

if (! -e $testvarfile) { 
	die "$testvarfile does not exist\n";
}
if (! -e $geneannotationlistvar) {
	die "$geneannotationlistvar does not exist\n";
}
if (! -e $dbnsfpfile) {
	die "$dbnsfpfile does not exist\n";
}



# load in gene annotations from file
print "Preloading gene annotations from $geneannotationlistvar\n";
if ($geneannotationlistvar =~ /.gz$/) {
	open (GENELISTVAR, "zcat $geneannotationlistvar |") or die "Cannot read $geneannotationlistvar: $!\n";
} else {
	open (GENELISTVAR, "$geneannotationlistvar") or die "Cannot read $geneannotationlistvar: $!\n";
}
<GENELISTVAR>;
while (<GENELISTVAR>) {
	$_ =~ s/\s+$//;
	loadGeneAnnot($_, \%geneannotations);
}
close GENELISTVAR;



# get header for dbNSFP annotations
# note that dbNSFP is 1-based not 0-based like CGI data
my $dbnsfp_headerline;
if ($dbnsfpfile =~ /.bz2$/) {
	$dbnsfp_headerline = `bzcat $dbnsfpfile | head -1`;
} else {
	$dbnsfp_headerline = `head -1 $dbnsfpfile.tsv`;
}	
$dbnsfp_headerline =~ s/\s+$//;
$dbnsfp_headerline =~ s/^#//;
my @dbnsfp_header = split("\t", $dbnsfp_headerline);
print "Reading dbNSFP from $dbnsfpfile\n";			
%dbnsfp_annots = ();							
$dbnsfp_annot_ref = loaddbNSFP($dbnsfp_annot_ref, \@dbnsfp_keepcol, $dbnsfpfile);



print "Annotating and calculating\n";		

open (OUT, ">$outputfile");
print OUT "chr\tbegin\tend\tvartype\tref\talt\txref\taltfreq\tMAF\tnhomref\tnhet\tnhomalt\tnmiss\tnhommaj\tnhommin\tCR\tgenesymbol\torientation\tcomponent\tcomponentIndex\thasCodingRegion\timpact\tnucleotidePos\tproteinPos\tannotationRefSequence\tsampleSequence\tgenomeRefSequence";
print OUT "\t".join("\t", @dbnsfp_header[@dbnsfp_keepcol])."\n";
if ($testvarfile =~ /.bz2$/) {
	open (FILE, "bzcat $testvarfile |") or die "Cannot read $testvarfile: $!\n";
} elsif ($testvarfile =~ /.gz/) {
	open (FILE, "zcat $testvarfile |") or die "Cannot read $testvarfile: $!\n";
} else {
	open (FILE, "$testvarfile") or die "Cannot read $testvarfile: $!\n";
}
my $header = <FILE>;
$header =~ s/\s+$//;
my @headerline = split("\t", $header);
my @ASMids = @headerline[8..$#headerline];
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($thisvarnum, $thischr, $thisstart, $thisend) = @line[0..3];
	$line[1] =~ s/chr//;														# make chromosome IDs just the number instead of "chr1", make it just "1" to make files more compatible with tabix
	my $vartype = $line[4];
	my ($ref, $alleleSeq) = @line[5..6];

	my $thiscoord = join('_', @line[1..6]);
	
	my @origgenotypes = @line[8..$#line];
	my @genotypes = @origgenotypes[0..2,4..76,78..$#origgenotypes];				# exclude 2 affected individuals (will throw off MAF calculations for Hutterite population frequency)
	
	my @genocounts = (0,0,0,0);
	for (my $subjidx=0; $subjidx<=$#genotypes; $subjidx++) {
		my $geno = $genotypes[$subjidx];
		if ($thischr =~ 'chrX') {
			if ($geno eq '00')  { $genocounts[0]++; }
			if ($geno eq '01' || $geno eq '10')  { $genocounts[1]++; }
			if ($geno eq '11')  { $genocounts[2]++; }			
			if ($geno =~ 'N')  { $genocounts[3]++; }
		} elsif ($thischr =~ 'chrM' || $thischr =~ 'chrY') {
			if ($geno eq '0')  { $genocounts[0]++; }
			if ($geno eq '1' || $geno eq '10')  { $genocounts[1]++; }
			if ($geno =~ 'N')  { $genocounts[3]++; }
		} else {
			if ($geno eq '00')  { $genocounts[0]++; }
			if ($geno eq '01' || $geno eq '10')  { $genocounts[1]++; }
			if ($geno eq '11')  { $genocounts[2]++; }
			if ($geno =~ 'N')  { $genocounts[3]++; }
		}
	}
	
	my $ngenosubj = $genocounts[0]+$genocounts[1]+$genocounts[2];	
	
	if ($ngenosubj != 0) {			# only print out lines that have at least one (high quality) genotyped individual!!
		my $altfreq = ($genocounts[1]+2*$genocounts[2])/(2*$ngenosubj);
		my $callrate = 1-($genocounts[3]/scalar(@genotypes));
	
		my ($maf, $genocounts_relhutt_maj, $genocounts_relhutt_min) = (0,0,0);		# count number of each genotype class using the hutterite minor allele as the alternative allele
		if ($genocounts[0] >= $genocounts[2]) {			# more or same ref alleles than alt alleles
			$maf = ($genocounts[1]+2*$genocounts[2])/(2*$ngenosubj);
			$genocounts_relhutt_maj = $genocounts[0];
			$genocounts_relhutt_min = $genocounts[2];
		} else {																	# same or more alt alleles as ref alleles
			$maf = ($genocounts[1]+2*$genocounts[0])/(2*$ngenosubj);
			$genocounts_relhutt_maj = $genocounts[2];
			$genocounts_relhutt_min = $genocounts[0];
		}

		my $geneannot = (("\t") x 10);
		if (exists $geneannotations{$thiscoord}) {
			$geneannot = $geneannotations{$thiscoord};
		} 
		
		# print output (combine genotypes with CGI gene annotations and dbNSFP annotations)
		# print OUT "$thisvarnum\t";
		for (my $i=1; $i<=7; $i++) {
			if ($line[$i] ne "") {
				print OUT "$line[$i]\t";
			} else {
				print OUT "\t";
			}
		}
		print OUT "$altfreq\t$maf\t".join("\t", @genocounts)."\t$genocounts_relhutt_maj\t$genocounts_relhutt_min\t$callrate\t$geneannot";
			
		my $dbnsfp_lookup = join("_", ($thischr, $thisend, $ref, $alleleSeq));
		if (($thisend-$thisstart)==1 && exists $dbnsfp_annots{$dbnsfp_lookup}) {		# haven't decided what to print out for ins/del/complex variants that span multiple bases		
			print OUT "\t$dbnsfp_annots{$dbnsfp_lookup}\n";		
		} else {
			print OUT "\t".(("\t") x scalar(@dbnsfp_keepcol))."\n";			
		}
		# END print output (combine genotypes with CGI gene annotations and dbNSFP annotations)
	} else {
		$countSNPnogenosubj++;
	}
}
close FILE;
close OUT;




sub loaddbNSFP {
	my $dbnsfp_annot_ref = $_[0];
	my $dbnsfp_keepcol_ref = $_[1];
	my $dbnsfpfile = $_[2];
	
	# prepare dbNSFP file for reading
	my $dbnsfp_handle;
	if (-e $dbnsfpfile) {
		open ($dbnsfp_handle, "bzcat $dbnsfpfile |") or die "Cannot read $dbnsfpfile: $!\n";
	} else {
		open ($dbnsfp_handle, "$dbnsfpfile") or die "Cannot read $dbnsfpfile: $!\n";
	}
	<$dbnsfp_handle>;				# skip header
	while (<$dbnsfp_handle>) {
		$_ =~ s/\s+$//;
		my @line = split("\t", $_);
		my ($thischr, $pos, $ref, $alt) = @line[0..3];
		my $dbnsfp_lookup = join("_", ("chr$thischr", $pos, $ref, $alt));
		# parse and save desired fields
		if (defined $line[0]) {
			${$dbnsfp_annot_ref}{$dbnsfp_lookup} = join("\t", @line[@{$dbnsfp_keepcol_ref}]);
		}
	}
	close $dbnsfp_handle;
	
	return $dbnsfp_annot_ref;
}

print "There were $countSNPnogenosubj SNPs where no subject had a high-quality genotype\n";



sub loadGeneAnnot {
	my ($string, $geneannot_ref) = @_;
	my @line = split ("\t", $_);
	my ($thischr, $thisstart, $thisend) = @line[0..2];

	my $lookup = join('_', @line[0..5]);
	
	# 0 chromosome	begin	end	varType	reference	call	xRef
	# 7 geneId	mrnaAcc	proteinAcc	symbol	orientation	component	componentIndex
	# 14 codingRegionKnown	impact	nucleotidePos	proteinPos	annotationRefSequence
	# 19 sampleSequence	genomeRefSequence
	
	my $annot = $line[10];
	for (my $i=11; $i<=20; $i++) {
		if ($line[$i]) {
			$annot .= "\t$line[$i]";
		} else {
			$annot .= "\t";
		}
	}
	${$geneannot_ref}{$lookup} = $annot;
}
                                                                                                                                                                           





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--testvar\ttestvariants file\n";
	print "\t--genelistvar\tgenelistvar file (produced by Generate_Gene_ListVariants_From_Gene_File_0_1_4.pl)\n";
	print "\t--dbnsfp\tfolder containing dbnsfp files, one per chr (from http://sites.google.com/site/jpopgen/dbNSFP)\n";
	print "\t--o\toutput <output.popinfo.tsv> file\n";
	die;
}