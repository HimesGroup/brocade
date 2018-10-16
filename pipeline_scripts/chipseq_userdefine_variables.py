#! /usr/bin/python

###
# Reference files
###

# hg38 reference files
hg38_fa = "/project/bhimeslab/Reference/hg38/genome.ERCC.fa"
hg38_bwa_index_prefix = "/project/bhimeslab/Reference/hg38/bwa_index/hg38"
hg38_len = "/project/bhimeslab/Reference/hg38/hg38.len"

# hg19 reference files
hg19_fa = "/project/bhimeslab/Reference/hg19/genome.ERCC.fa"
hg19_bwa_index_prefix = "/project/bhimeslab/Reference/hg19/bwa_index/hg19"
hg19_len = "/project/bhimeslab/Reference/hg19/hg19.len"

# mm38 reference files
mm38_fa = "/project/bhimeslab/Reference/mm38/mm.GRCm38.genome.ERCC92.fa"
mm38_bwa_index_prefix = "/project/bhimeslab/Reference/mm38/bwa_index/mm38"
mm38_len = "/project/bhimeslab/Reference/mm38/mm38.len"

# mm10 reference files
mm10_fa = "/project/bhimeslab/Reference/mm38/UCSC_mm10_genome.ERCC92.fa"
mm10_bwa_index_prefix = "/project/bhimeslab/Reference/mm10/bwa_index/mm10"
mm10_len = "/project/bhimeslab/Reference/mm10/mm10.len"

# rn6 reference files
rn6_fa = "/project/bhimeslab/Reference/rn6/genome.ERCC.fa"
rn6_bwa_index_prefix = "/project/bhimeslab/Reference/rn6/bwa_index/rn6"
rn6_len ="/project/bhimeslab/Reference/rn6/rn6.len"

# susScr3 reference files
susScr3_fa = "/project/bhimeslab/Reference/susScr3/Ensembl_susScr3_genome.ERCC92.fa"
susScr3_bwa_index_prefix = "/project/bhimeslab/Reference/susScr3/bwa_index/susScr3"
susScr3_len = "/project/bhimeslab/Reference/susScr3/susScr3.len"

###
# BigData UTL
###

# user-provided public URL to put bigwig files
bigdata_url="http://public.himeslab.org"

###
# Informatics Tools
###

# path and directory
trimmomatic="/opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar"
picard_dir="/opt/software/picard/picard-tools-1.96/"

# version
bwa_version="0.7.10-r789"
trimmomatic_version="0.32"
fastqc_version="0.11.7"
samtools_version="1.8"
bamtools_version="2.3.0"
picard_version="1.96"
macs2_version="2.1.1"
deseq2_version="1.18.1"
diffbind_version="2.6.6"
chipseeker_version="1.14.2"

###
# Genes of interest
###

fav_gene=["CRISPLD2","ORMDL3"]

###
# Other information
###

# author name and contact info (shown in .html report)
author="Mengyuan Kan (mengykan@upenn.edu)"
