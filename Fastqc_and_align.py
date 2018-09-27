#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 11:21:55 2018

@author: Avantika Diwadkar (avantika.diwadkar@pennmedicine.upenn.edu)
"""

#Libraries 
import os
import subprocess
import urllib
import urllib.request
from SRAdownload import *


#Set up Reference genome files (if not already set up)
ref_path = = input("Please enter reference genome full path: ")

#Perform FASTQC on the sample files
fastq_cmd = "fastqc %s/fastq/*.fastq.gz  -o %s/"%(geo_dir,res_dir)

#Alignment of reads to genome
bwa_dir = data_dir + "/bwa"
bam_dir = data_dir + "/bam"
fastq_dir = geo_dir + "/fastq"


Bcall = '''
fastq_dir="%s"
bam_dir="%s"

# use samples for GR_dex vs GR comparison
ls $fastq_dir | sed 's/.fastq.gz//g' > samples.txt
threat=12
while read samp; do
echo "bwa mem -t $threat %s/hg38.fa $fastq_dir/$samp.fastq.gz | samtools view -S -b - | samtools sort -@$threat -T $bam_dir/$samp.tmp -o $bam_dir/$samp.bam -" # alignment and sort$
echo "samtools index $bam_dir/$samp.bam" # index
done < samples.txt | mk_bsub.py --line 2 --thread 12 --job_name chipseq_bwa
'''%(fastq_dir,bam_dir,ref_path)    


if __name__ == '__main__':

	#Fastqc on sample files
	print("Fastqc in progress")
	#subprocess.call(fastq_cmd,shell=True)

	#Make bam and bwa folders in data directory
	make_dir(bam_dir)    
	make_dir(bwa_dir)

	#Write alignment file
	write_in_file("alignment.sh",Bcall)

	run_cmd = "chmod +x alignment.sh; ./alignment.sh"
	subprocess.call(run_cmd, shell=True)
	print("Alignment file run complete - lsf files created")


