#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import glob

# import user-defined parameters
import chipseq_userdefine_variables as userdef # read in user-defined variable python script
# import HPC parameters
memory=userdef.memory
queue=userdef.queue

def check_exist(file):
    if not os.path.exists(file):
        print "The file "+file+" does not exist. Please check!"
        sys.exit()

def get_sample_info(fin):
    """
    Read in information from phenotype file
    Create a dictionary to store each column
    """
    f = open(fin,'r')
    f = f.read().split('\n')
    f = map(lambda x: x.rstrip(), f)

    if '' in f:
        f.remove('')
    header = f[0].split('\t') # list
    c = f[1:]

    # Obtain column index from file header
    if "Sample" not in header:
        print "Sample column is missing. Please check!"
        sys.exit()

    # Obtain column index from file header
    if "Antibody" not in header:
        print "Antibody column is missing. Please check!"
        sys.exit()

    d = {} # key: column name, value: column

    for i in range(len(header)):
        colname=header[i]
        d[colname]=map(lambda x: x.split('\t')[i],c)

    return d

def macs2_peakcall(curr_sample, curr_control, curr_peak, path_start):
    """
    Use macs2 for peak calling
    """

    control_bam=path_start+curr_control+"/bwa_out/"+curr_control+".bam"
    treatment_bam=path_start+curr_sample+"/bwa_out/"+curr_sample+".bam"
    macs2_dir=path_start+"/"+curr_sample+"/"+"macs2_out/"
    cmd="mkdir "+macs2_dir+"\n"
    if curr_peak=="narrow":
        cmd=cmd+"macs2 callpeak -c "+control_bam+" -t "+treatment_bam+" -n "+curr_sample+" --outdir "+macs2_dir+" -f BAM -g hs -B -q 0.01"
    else:
        cmd=cmd+"macs2 callpeak -c "+control_bam+" -t "+treatment_bam+" -n "+curr_sample+" --outdir "+macs2_dir+" --broad -f BAM -g hs -B -q 0.01"

    cmd=cmd+"\n"
    return cmd


def blacklist_rm(curr_sample, curr_peak, ref_genome, path_start, template_dir):
    """
    Remove peaks within blacklisted regions
    """

    macs2_dir=path_start+curr_sample+"/macs2_out/"

    if curr_peak=="narrow":
        peak_fn=macs2_dir+curr_sample+"_peaks.narrowPeak"
    else:
        peak_fn=macs2_dir+curr_sample+"_peaks.broadPeak"

    if ref_genome=="hg38" or ref_genome=="hg19":

        # check if blacklisted region file exists
        blacklist_fn = template_dir+ref_genome+"_blacklist.bed"
        if not os.path.exists(blacklist_fn):
            print "Cannot find blacklisted region file: "+blacklist_fn
            print "Please download it from https://github.com/HimesGroup/brocade"
            sys.exit()

        cmd="bedtools intersect -v -a "+peak_fn+" -b "+blacklist_fn+" | grep -P 'chr[\dXY]+[ \\t]' | awk 'BEGIN{OFS=\"\\t\"} {if ($5>1000) $5=1000; print $0}' > "+macs2_dir+curr_sample+".blackfilt.bed"
    else: # no black list for other genomes
        cmd="cat "+peak_fn+" | grep -P 'chr[\dXY]+[ \\t]' | awk 'BEGIN{OFS=\"\\t\"} {if ($5>1000) $5=1000; print $0}' > "+macs2_dir+curr_sample+".blackfilt.bed"

    cmd=cmd+"\n"
    return cmd


def merge_input_bam(info_dict,path_start):
    """
    Find the input files in the sample info sheet and merge them into one single bam file
    """
    input_files = []
    NA_index = [i for i, x in enumerate(info_dict["Input"]) if x == "NA"]
    all_files = info_dict["Sample"] #full file name with path
    for i in NA_index:
        file = path_start+all_files[i]+"/bwa_out/"+all_files[i]+".bam"
        input_files.append(file)

    ffiles = " ".join(str(x) for x in input_files)
    return ffiles


def lsf_file(job_name, cmd, memory=memory, thread=1, queue=queue):
    """
    Creates .lsf files
    """

    outp = open(job_name+".lsf",'w')
    outp.write("#!/bin/bash\n")
    outp.write("#BSUB -L /bin/bash\n")
    outp.write("#BSUB -J "+job_name+"\n")
    outp.write("#BSUB -q "+queue+"\n")
    outp.write("#BSUB -o "+job_name+"_%J.out\n")
    outp.write("#BSUB -e "+job_name+"_%J.screen\n")
    outp.write("#BSUB -M "+str(memory)+"\n")
    outp.write("#BSUB -n "+str(thread)+"\n")
    outp.write(cmd)
    outp.write("\n")
    outp.close()


def main(sample_info_file, project_name, ref_genome, path_start, template_dir,input_bam):
    """
    Peak calling using macs2
    Remove peaks within blacklisted regions
    """

    ####
    # Set up and check
    ####

    # Set up project directory
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
	path_start = path_start+"/"

    # Set up template directory
    if template_dir == "./":
        template_dir = os.getcwd()
    if template_dir[-1] != "/":
	template_dir = template_dir+"/"


    # check if sample info exists
    if not os.path.exists(sample_info_file):
        print "Cannot find sample_info_file: "+sample_info_file
        sys.exit()

    # Obtain sample information from phenotype file
    print "sample_info_file = " + sample_info_file
    check_exist(sample_info_file)
    info_dict = get_sample_info(sample_info_file)
    sample_names = info_dict["Sample"]

    # Obtain input DNA sample
    if "Input" not in info_dict:
        print "No input DNA is provided."
        sys.exit()

    controls = info_dict["Input"]

    #Concatenate the input files if input_bam
    if input_bam:
        print "Merging input DNA control .bam files. It will take a while..."
        ffiles = merge_input_bam(info_dict,path_start)
        input_path = path_start+project_name+'_input'
        mk_dir = 'mkdir -p '+input_path+'/bwa_out/'
        input_merge = 'samtools merge ' + input_path+'/bwa_out/'+project_name+"_input.bam "+ ffiles
        subprocess.Popen(mk_dir,shell=True).wait()
        subprocess.Popen(input_merge,shell=True).wait()
        print "Done!"

    # remove NA
    control_idx=filter(lambda x: controls[x] not in 'NA', range(len(sample_names)))

    # Obtain peak information
    if "Peak" not in info_dict:
        print "Peak column that specifies narrow or broad peak is missing. Please check!"

    peaks = info_dict["Peak"]
    peak_filt=map(lambda x: peaks[x], control_idx)
    if any(map(lambda x:x not in ['narrow','broad'], peak_filt)):
        print "Specify peaks as 'narrow' or 'broad'. Please check!"
        sys.exit()

    for i in control_idx:

        cmd=""
        curr_sample=sample_names[i]
        curr_control=controls[i]
        curr_peak=peaks[i]
        print "treatment and control pair: "+curr_sample+"+"+curr_control

        ###
        # Peak calling
        ###

        macs2_cmd=macs2_peakcall(curr_sample, curr_control, curr_peak, path_start)
        cmd=cmd+macs2_cmd

        ###
        # Remove blacklisted regions
        ###

        blackrm_cmd=blacklist_rm(curr_sample, curr_peak, ref_genome, path_start, template_dir)
        cmd=cmd+blackrm_cmd

        ###
        # Create .lsf files
        ###
        lsf_file(curr_sample+"_macs2", cmd, thread=1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a project.")
    parser.add_argument("--project_name", type=str, help="Prefix name of for output directory and files")
    parser.add_argument("--samples_in", help="A tab-delimited txt file containing sample information with full path. See example file: sample_info_file.txt")
    parser.add_argument("--ref_genome", default="hg38", type=str, help="Specify reference genome to filter out blacklisted regions (options: hg38, hg19, mm38, mm10, rn6)")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
    parser.add_argument("--template_dir", default="./", type=str, help="directory to put provided or user defined reference index files")
    parser.add_argument("--input_bam", action='store_true', help="If specified, merges all input bam files into one file.")
    args = parser.parse_args()

    if args.project_name is None or args.samples_in is None:
        parser.print_help()
        sys.exit()

    main(args.samples_in, args.project_name, args.ref_genome, args.path_start, args.template_dir,args.input_bam)

