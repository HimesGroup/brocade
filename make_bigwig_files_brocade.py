#!/usr/bin/python
import sys
import subprocess
import os
import math
import string
import argparse


def make_track_file(path_start,project):
    #Make text file with all bdgs file names 
    bdg_cmd = "cd %s/%s/data/macs2/bedGraph ; ls *.bdg > %s/sample_bdgs.txt ;cd %s "%(path_start,project,path_start,path_start)
    subprocess.call(bdg_cmd, shell=True)
    with open(path_start+"/sample_bdgs.txt") as f:
        files = f.read().splitlines()
     
    #Get conditions of diffbinding
    conditions = []
    for i in files:
         Ab = i.split("_")
         conditions.append(Ab[0])     
     
    conditions = list(set(conditions))
     
    #get list of labels
    labels=[]
    for i in files:
        lab = i.split(".bdg")
        labels.append(lab[0])
    	
    	#Get a dictionary of colors for conditions in info_sheet
    colors = ["0,0,400", "400,0,0", "0,0,100", "100,0,0", "200,200,0", "0,200,200"]
    color_dic = {}
    for i in range(len(conditions)):
    		color_dic[conditions[i]] = colors[i]
      
    #Make track file and set url link
    outp = open(path_start+"/ucsc_track_names.txt", "w")
    output_url = "http://public.himeslab.org/"
    
    
    	#Write out each line of the track names file
    count=0
    for k in files:
        curr_bigwig = output_url+project+"/"+labels[count]+".bw"
        outp.write("track type=bigWig name=\"")
        outp.write(labels[count]+"\" ")
        cond = (labels[count].split("_"))[0]
        outp.write("color="+color_dic[cond])
        #outp.write(" gridDefault=on maxHeightPixels=50 visibility=full autoScale=off viewLimits=0:13000 description=")
        outp.write(" gridDefault=on maxHeightPixels=50 visibility=full autoScale=on description=")
        outp.write(labels[count]+"\" bigDataUrl=")
        outp.write(curr_bigwig+"\n")
        count += 1
    outp.close()
    print ("Created text file to load tracks in UCSC genome browser")
    return(files,labels)





def main(project_name,ref_genome):
    """
	Creates sorted bedgraph and bigwig files of all samples for a batch for visualization in UCSC genome browser
	"""
    path_start = os.getcwd()
    if path_start[-1] != "/":
        path_start = path_start+"/"
        
    out_dir = path_start+project_name+"/"+"data/macs2/"

    bedgraph_out = out_dir+"bedGraph/"

    bigwig_out = bedgraph_out+"bigwig/"
    #if not os.path.exists(bigwig_out):
        #os.makedirs(bigwig_out)
		
    runs= make_track_file(path_start,project_name)
    files, labels = runs

    for k in range(len(labels)):
        job_name = labels[k]
        curr_sample = files[k]
        
		
        if ref_genome == "hg38":
            chr_size_file = "/project/bhimeslab/Reference/hg38/hg38.chrom.sizes"
        elif ref_genome == "mm38":
            chr_size_file = "/project/bhimeslab/Reference/mm38/mm38.chrom.sizes"		
        elif ref_genome == "mm10":
            chr_size_file = "/project/bhimeslab/Reference/mm38/mm10.chrom.sizes"	
        elif ref_genome == "rn6":
            chr_size_file = "/project/bhimeslab/Reference/rn6/rn6.chrom.sizes"			
        elif ref_genome == "susScr3":
            chr_size_file = "/project/bhimeslab/Reference/susScr3/susScr3.chrom.sizes"	
        else:
            print ("No chromosome size file for this genome available:", ref_genome)
		
		#Make lsf file		
        outp = open(job_name+".lsf", "w")
        outp.write("#!/bin/bash \n")
        outp.write("#BSUB -L /bin/bash\n")
        outp.write("#BSUB -J "+job_name+"\n")
        outp.write("#BSUB -q normal \n")
        outp.write("#BSUB -o "+job_name+"_%J.out\n")
        outp.write("#BSUB -e "+job_name+"_%J.screen\n")
        outp.write("#BSUB -M 36000\n")
        outp.write("#BSUB -n 12\n")
		
		#Sort bedgraph file
        outp.write("LC_COLLATE=C sort -k1,1 -k2,2n "+bedgraph_out+curr_sample+" > "+bedgraph_out+job_name+".sorted.bdg\n")

		#Create bigwig files
        outp.write("bedGraphToBigWig "+bedgraph_out+job_name+".sorted.bdg "+chr_size_file+" "+bigwig_out+job_name+".bw\n")

        outp.close()
	   #subprocess.call("bsub < "+job_name+".lsf", shell=True)
	   #subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creates sorted bedgraph and bigwig files for ChIP-seq samples associated with a project.")
    #parser.add_argument("--path_start", default="./", type=str, help="Directory path to where project directory resides (default=./)")
    	#Change this to include chr size files into main directories with ref files (see main too)
    	#parser.add_argument("--chr_size_file", default="/project/bhimeslab/Reference/hg38/hg38.chrom.sizes", type=str, help="A file created using UCSC browser script fetchChromSizes (e.g. fetchChromSizes hg38 > hg38.chrom.sizes). Default: /project/bhimeslab/Reference/hg38/hg38_ERCC.chrom.sizes")
    #parser.add_argument("--aligner", default="", type=str, help="Should TopHat, STAR, or bowtie2 be used as aligner?"
    		#"(options: tophat, star, bowtie2)")
    parser.add_argument("--project_name", type=str, help="Name of project that all samples correspond to.")
    parser.add_argument("--ref_genome", type=str, help="Reference genome used for peak calling.")
    #parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
    args = parser.parse_args()
    main(args.project_name, args.ref_genome)



