#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 11:21:55 2018

@author: diwadkar
"""

#Libraries 
import os
import subprocess
import urllib
import urllib.request
from SRAdownload import *



#Perform FASTQC on the sample files
fastq_cmd = "fastqc %s/fastq/*.fastq.gz  -o %s/"%(geo_dir,res_dir)


#Get human reference genome files -ERROR
ref_path = '/project/bhimeslab/Reference/hg38/' 
hg_url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
g_cmd = "gunzip hg38.fa.gz"


#Build human genome
bwa_cmd = "bwa index -a %sbwtsw hg38.fa -p hg38"%(ref_path)

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
echo "bwa mem -t $threat /project/bhimeslab/Reference/hg38/hg38.fa $fastq_dir/$samp.fastq.gz | samtools view -S -b - | samtools sort -@$threat -T $bam_dir/$samp.tmp -o $bam_dir/$samp.bam -" # alignment and sort$
echo "samtools index $bam_dir/$samp.bam" # index
done < samples.txt | mk_bsub.py --line 2 --thread 12 --job_name chipseq_bwa
'''%(fastq_dir,bam_dir)    


if __name__ == '__main__':

	#Fastqc on sample files
	print("Fastqc in progress")
	#subprocess.call(fastq_cmd,shell=True)

	#Set up hg38 if not already present
	#download_file(hg_url, ref_path)
	print("Downloaded hg38.fa.gz")
	#subprocess.call(g_cmd,shell=True)
	print("Unzipped the file")
	#subprocess.call(bwa_cmd, shell=True)
	print("Aligning to genome hg38")


	#Make bam and bwa folders in data directory
	make_dir(bam_dir)    
	make_dir(bwa_dir)


	#Write alignment file
	write_in_file("alignment.sh",Bcall)

	run_cmd = "chmod +x alignment.sh; ./alignment.sh"
	#subprocess.call(run_cmd, shell=True)
	print("Alignment file run complete - lsf files created")

	lsf_cmd = "for i in chipseq_bwa*.lsf; do echo $i; bsub < $i; done"
	#subprocess.call(lsf_cmd, shell=True)
	print("Jobs submitted")

