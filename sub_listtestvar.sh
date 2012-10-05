#!/bin/bash
#$ -cwd


# use listvariants to make file from all masterVar files listed in text file

~/bin/cgatools listvariants --beta --reference /ober_resources/obertools/lib/CGAtools/build37.crr --variants `cat /clusta/jxchong/CGI_WGS/all.mastervarpath.2012-09-20.txt` --output /clustb/jxchong/CGI_WGS/all.2012-09-20.listvar.tsv 

cp /clustb/jxchong/CGI_WGS/all.2012-09-20.listvar.tsv /clusta/jxchong/CGI_WGS/all.2012-09-20.listvar.tsv

# use testvariants to make testvar file

~/bin/cgatools testvariants --beta --reference /ober_resources/obertools/lib/CGAtools/build37.crr --input /clusta/jxchong/CGI_WGS/all.2012-09-20.listvar.tsv --variants `cat /clusta/jxchong/CGI_WGS/all.mastervarpath.2012-09-20.txt` --output /clustb/jxchong/CGI_WGS/all.2012-09-20.testvar.tsv

cp /clustb/jxchong/CGI_WGS/all.2012-09-20.testvar.tsv /clusta/jxchong/CGI_WGS/all.2012-09-20.testvar.tsv
