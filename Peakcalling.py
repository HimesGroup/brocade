#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 11:48:34 2018

@author: diwadkar
"""

#Libraries 
import os
import subprocess
import urllib
import urllib.request
from SRAdownload import *
from Fastqc_and_align import *



#Add SRR sample IDs in the Phenotype file
sample_file = "/home/diwadkar/samples.txt"

f=open(sample_file,"r")
flines=f.readlines()
#flines = lines[1:]
result=[]
for x in flines:
    result.append(x.split('\n')[0])
f.close()


phenotype_file = res_dir+"/"+geo_id+"_Phenotype_withoutQC.txt"
new_phenotype_file = res_dir+"/"+geo_id+"_Phenotype_withoutQC_withSRR.txt"

#Peak calling using MACS2
peak_dir = data_dir + "/macs2"


Pcall = '''
bam_dir="%s"
macs2_dir="%s"

echo "macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309351.bam $bam_dir/SRR5309353.bam -n GR_control --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #GR_Control
macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309351.bam -n GR_1D --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309353.bam -n GR_1E --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309354.bam $bam_dir/SRR5309356.bam -n GR_dex --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #GR_Dex
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309354.bam -n GR_2A --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309356.bam -n GR_2E --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309357.bam $bam_dir/SRR5309358.bam -n GR_dex_GFP --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #GR_dex_GFP
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309357.bam -n GR_3D --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309358.bam -n GR_3E --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309359.bam $bam_dir/SRR5309360.bam -n GR_dex_KLF15 --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #GR_dex_KLF15
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309359.bam -n GR_4B --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309360.bam -n GR_4C --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309361.bam $bam_dir/SRR5309362.bam -n RNAP2_control --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #RNAP2_Control
macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309361.bam -n RNAP2_5E --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309362.bam -n RNAP2_5b --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309363.bam $bam_dir/SRR5309364.bam -n RNAP2_Dex --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #RNAP2_Dex
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309363.bam -n RNAP2_6d --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309364.bam -n RNAP2_6e --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309365.bam $bam_dir/SRR5309366.bam -n RNAP2_dex_GFP --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #RNAP2_dex_GFP
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309365.bam -n RNAP2_7c --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309366.bam -n RNAP2_7e --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309367.bam $bam_dir/SRR5309368.bam -n RNAP2_control_KLF15 --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #RNAP2_control_KLF15
macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309367.bam -n RNAP2_8B --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309352.bam -t $bam_dir/SRR5309368.bam -n RNAP2_8D --outdir $macs2_dir -f BAM -g hs -B -q 0.01

macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309369.bam $bam_dir/SRR5309370.bam -n RNAP2_dex_KLF15 --outdir $macs2_dir -f BAM -g hs -B -q 0.01 #RNAP2_dex_KLF15
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309369.bam -n RNAP2_9b --outdir $macs2_dir -f BAM -g hs -B -q 0.01
macs2 callpeak -c $bam_dir/SRR5309355.bam -t $bam_dir/SRR5309370.bam -n RNAP2_9e --outdir $macs2_dir -f BAM -g hs -B -q 0.01" | /home/diwadkar/mk_bsub.py --line 1 --thread 1 --job_name macs2 
'''%(bam_dir,peak_dir)


if __name__ == '__main__':

	#Add SRR numbers to phenotype file
	#if OS X then add '' -e after -i in sed
	srr_cmd = "sed -i $'1 i\\\nSRR' %s ; paste %s %s > %s"%(sample_file,sample_file,phenotype_file,new_phenotype_file)
	subprocess.call(srr_cmd, shell=True)

	#Make peak calling directory
	make_dir(peak_dir)

	#Make peakcalling file and run
	write_in_file("peakcalling.sh",Pcall)
	subprocess.call("chmod +x peakcalling.sh; ./peakcalling.sh", shell=True)   
	#subprocess.call("for i in macs2*.lsf; do echo $i; bsub < $i; done", shell=True)

	#Convert narrowPeak to bed files 
	bedcall = "cd %s; for i in *.narrowPeak; do echo $i ; cut -f 1-6 $i > $i.bed; done"%(peak_dir)
	#subprocess.call(bedcall, shell=True)


