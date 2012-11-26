General information:
*   popinfo file contains population info such as: MAF; alternate allele freq; genotype counts; call rate; gene-based annotations via CGI (if it's in a gene, what the gene name is, what part of the gene; if it's a coding variant, whether it's missense, nonsense, stopgain, stoploss, splicing; what exon or intron it's in; which # amino acid is affected
*   all numbers based on genomic coordinates (such as the exon number affected by a given variant, or the variant's position) are 0-based not 1-based, i.e. counting starts from exon 0
*   popinfo file integrates annotations such as SIFT and GERP scores from dbNSFP but dbNSFP is not necessary.  dbNSFP (the most recent iteration I used was 2.0b3) seemed to be missing scores for some SNVs that ought to have SIFT/GERP scores; not sure wh

To generate popinfo file:

1.   sub_listtestvar.sh (example submission script)
	1.  Use cgatools listvariants to generate file with all variants in the genomes in question (I limited it to VQHIGH variants only but this is optional).
 	2.  Use cgatools testvariants to generate file with "genotypes" (test for presence of a given variant in each genome).
2.   sub_genelistvar.sh (example submission script)
	1.  Use CGI's unofficial Generate_Gene_ListVariants_From_Gene_File_0_1_4 to get gene-based annotations for variants.  Unfortunately cgatools testvariants and listvariants don't retain the gene-based annotation information!!!
3.   Run split_by_chr.pl to split the files generated in steps 1 and 2 by chromosome (for faster processing/batching)
4.   runsub_popinfo.sh (example submission script): Generate popinfo file
5.   Not implemented yet, but if samtools-associated programs bgzip/tabix is installed, should use that to compress the popinfo file for quick/better lookup (instead of reading through the popinfo file line by line).
