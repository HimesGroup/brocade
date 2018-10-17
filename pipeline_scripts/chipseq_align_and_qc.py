#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import glob


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

    d = {} # key: column name, value: column

    for i in range(len(header)):
        colname=header[i]
        d[colname]=map(lambda x: x.split('\t')[i],c)

    return d



def get_genome_ref_files(genome):
    """
    Location of all reference files needed for a given genome.
    Current choice: "hg38"
    """

    import chipseq_userdefine_variables as userdef  # read in user-defined variable python script: improt reference files

    if genome == "hg38":
	fa = userdef.hg38_fa
	bwa_index_prefix = userdef.hg38_bwa_index_prefix

    elif genome == "hg19":
	fa = userdef.hg19_fa
	bwa_index_prefix = userdef.hg19_bwa_index_prefix

    elif genome == "mm38":
	fa = userdef.mm38_fa
	bwa_index_prefix = userdef.mm38_bwa_index_prefix

    elif genome == "mm10":
	fa = userdef.mm10_fa
	bwa_index_prefix = userdef.mm10_bwa_index_prefix

    elif genome == "rn6":
	fa = userdef.rn6_fa
	bwa_index_prefix = userdef.rn6_bwa_index_prefix

    elif genome == "susScr3":
	fa = userdef.susScr3_fa
	bwa_index_prefix = userdef.susScr3_bwa_index_prefix

    else:
	print 'Unknown genome selected: ', genome
	sys.exit()

    # create reference bwa index files
    bwa_index_files=map(lambda x: bwa_index_prefix+x, ['.amb','.ann','.bwt','.pac','.sa'])

    # check if reference files exist
    map(lambda x:check_exist, [fa, bwa_index_files])

    return(fa, bwa_index_prefix)

def get_adapter_info(fin, index_type):
    """
    Read in information from adapter sequence file
    Obtain adapter sequences of corresponding index type
    """
    f = open(fin,'r')
    f = f.read().split('\n')
    f = map(lambda x: x.rstrip(), f)

    if '' in f:
        f.remove('')

    header = f[0].split('\t') # list
    c = f[1:]

    # check if Type, Index, Description, and Sequence are in the file
    if set(["Type", "Index", "Description", "Sequence"]) != set(header):
        print "Not all the following colunms Type, Index, Description, and Sequence are in the file. Please check!"
        sys.exit()

    d = {} # key: column name, value: column

    for i in range(len(header)):
        colname=header[i]
        d[colname]=map(lambda x: x.split('\t')[i],c)

    # check if index_type within provided types
    if index_type not in d["Type"]:
        print index_type+" is not in the provided index types: "+', '.join(map(str,list(set(d["Type"]))))
        print "Please check!"
        sys.exit()

    # obtain idx of corresponding index type
    idx=filter(lambda x:index_type in d["Type"][x], range(len(c)))

    # Create dictionary to store index sequences of corresponding index type
    d_seq = {} # key: Index; value: Sequence

    for i in idx:
        index_type=d["Type"][i]
        index_seq=d["Index"][i]
        index_des=">"+d["Description"][i] # add '>' before description
        index_allseq=d["Sequence"][i]

        if index_seq not in d_seq:
            d_seq[index_seq]=index_des+"\n"+index_allseq
        else:
            d_seq[index_seq]=d_seq[index_seq]+"\n"+index_des+"\n"+index_allseq

    return d_seq

def index_check(indexes, index_type, template_dir):
    """
    For user specified index, check whether they are in the provided adapter sequence files.
    """	

    indexes=[x for x in indexes if x !='NA'] # remove NA values

    if index_type is None:
        print "Index type (--index_type) is not specified. Please check."
        sys.exit()

    # if unique dual (UD) index adapters are used, the two indexes i7 and i5 combined by "+" are needed
    if index_type in ["illumina_ud_sys1", "illumina_ud_sys1"]:
        if not all(map(lambda x:'+' in x, indexes)):
            print "Unique dual (UD) index type is specified. Provide i7 and i5 sequences using i7+i5 in Index column. Please check!"
            sys.exit()

    index_fn=template_dir+"adapter_primer_sequences.txt"
    print "index_type = "+index_type
    print "Use provided adapter and primer sequence file: "+index_fn
    check_exist(index_fn)

    # obtain adapter and primer sequences based on index type and index sequence
    index_dict=get_adapter_info(index_fn, index_type)
    index_seqs=index_dict.keys()

    # check if index sequences within specified index type
    # check if all indexes sequences from phenotype file are in provided reference file
    index_diff=list(set(indexes)-set(index_seqs))
    if len(index_diff) >0:
        print "The index(es) from phenotype file are not in the provided adapter sequence file: "+', '.join(map(str,index_diff))
        print "Please check!"
        sys.exit()

    return index_dict

def make_adapter_fa(curr_sample, out_dir, curr_index, index_dict):
    """
    Make a .fa file with adapter and primer sequences of corresponding sample
    """
    fa_out=index_dict[curr_index] # output index description and sequences from index dictionary
    primer_keys=filter(lambda x:"Primer" in x, index_dict.keys()) # obtain Read1 and/or Read2 primer names from dictionary
    for key in primer_keys:
        fa_out=fa_out+"\n"+index_dict[key]
    outp=open(out_dir+curr_sample+"_adapter.fa",'w') # output in curr_sample+"_adapter.fa"
    outp.write(fa_out)
    outp.write("\n")
    outp.close()

def trim_and_fastqc(curr_sample, curr_index, out_dir, R1, R2, R1_trim, R2_trim):
    """
    Use Trimmomatic to trim adaptor and perform fastqc
    """
    
    # Trim adapter
    import chipseq_userdefine_variables as userdef # read in user-defined variable python script: import trimmomatic java path
    trimmomatic=userdef.trimmomatic

    cmd="" # create command variable
    if curr_index == 'NA':
        print curr_sample+": no index defined. Skip trimming."
    else:
        if R2=="":
            if os.path.isfile(R1_trim):
                print curr_sample+" trimmed file already exists. Skip trimming."
            else:
                cmd = "java -Xmx1024m  -classpath "+trimmomatic+" org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n"

        else:
            R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"
            if os.path.isfile(R1_trim) and os.path.isfile(R1_trim):
                print curr_sample+" R1 and R2 trimmed files already exist. Skip trimming."
            else:
                cmd = "java -Xmx1024m  -classpath " +trimmomatic+" org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n"

    # Fastqc after trimming
    # create standard fastqc .zip file name
    R1_fastqc_fn = out_dir+curr_sample+"_R1_Trimmed_fastqc.zip"
    R2_fastqc_fn = out_dir+curr_sample+"_R2_Trimmed_fastqc.zip"

    if (R2=="" and os.path.isfile(R1_fastqc_fn)) or (os.path.isfile(R1_fastqc_fn) and os.path.isfile(R2_fastqc_fn)):
        print curr_sample+" fastqc results already exist. Skip fastqc."
    else:
        # Run fastqc
        cmd=cmd+"fastqc -o "+out_dir+" "+R1_trim+" "+R2_trim+"\n" # since R2 could be "", no harm to add an empty string directly here
        if curr_index=='NA': # use no-trimmed files
            # original fastqc .zip file name
            R1_org_name=R1_trim.split(".fastq", 1)[0]+"_fastqc.zip"
            cmd=cmd+"cp "+R1_org_name+" "+R1_fastqc_fn+"\n"
            if R2_trim!="": # if R2 exists
                R2_org_name=R2_trim.split(".fastq", 1)[0]+"_fastqc.zip"
                cmd=cmd+"cp "+R2_org_name+" "+R2_fastqc_fn+"\n"

    #Get total number of reads, unique reads, % unique reads from trimmed file(s).
    if ".gz" in R1_trim: # for gzip .fastq file
        cmd=cmd+"zcat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+out_dir+curr_sample+"_ReadCount\n"
    else:
        cmd=cmd+"cat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+out_dir+curr_sample+"_ReadCount\n"
    if R2_trim!="":
        if ".gz" in R2_trim:
            cmd=cmd+"zcat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+out_dir+curr_sample+"_ReadCount\n"
        else:
            cmd=cmd+"cat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+out_dir+curr_sample+"_ReadCount\n"

    return cmd

def bwa_and_bamstat(curr_sample, out_dir, library_type, R1_trim, R2_trim, bwa_index_prefix):
    """
    Use BWA for alignment
    """

    cmd="mkdir "+out_dir+"bwa_out\n"
    cmd=cmd+"cd "+out_dir+"bwa_out\n"
    if R2_trim!="": # paired-end
        cmd=cmd+"bwa mem -t 12 "+bwa_index_prefix+" "+R1_trim+" "+R2_trim+" | samtools view -S -b - | samtools sort -@12 -T "+curr_sample+".tmp -o "+curr_sample+".bam -"
    else: # single-end
        cmd=cmd+"bwa mem -t 12 "+bwa_index_prefix+" "+R1_trim+" | samtools view -S -b - | samtools sort -@12 -T "+curr_sample+".tmp -o "+curr_sample+".bam -"
    cmd=cmd+"\n"

    #Create indexed bam file:
    cmd=cmd+"samtools index -@12 "+curr_sample+".bam\n"

    """
    Obtain QC metrics from bam file
    """
    # bam stat
    import chipseq_userdefine_variables as userdef # read in user-defined variable python script: import picard directory
    picard_dir=userdef.picard_dir

    #Write out index stats of where reads align to by chr:
    cmd=cmd+"samtools idxstats "+curr_sample+".bam > "+curr_sample+".stats\n"
    #Write out bamtools summary stats:
    cmd=cmd+"bamtools stats -in "+curr_sample+".bam > "+curr_sample+".bamstats\n"

    #Gather metrics unique to paired-end samples using CollectInsertSizeMetrics
    if library_type in ["PE"]:
	cmd=cmd+"java -Xmx2g -jar "+picard_dir+"CollectInsertSizeMetrics.jar VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE="+curr_sample+".InsertSizeHist.pdf INPUT="+curr_sample+".bam OUTPUT="+curr_sample+".InsertSizeMetrics\n"

    return cmd

def bam2bw_and_track(curr_sample, curr_color, out_dir, template_dir, track_fn, bigdata_path, len_fn):
    """
    Convert bam to bw file
    """

    cmd=""
    cmd=cmd+"genomeCoverageBed -split -bg -ibam "+curr_sample+".bam -g "+len_fn+" > "+curr_sample+".bdg\n"
    cmd=cmd+"LC_COLLATE=C sort -k1,1 -k2,2n "+curr_sample+".bdg > "+curr_sample+".sorted.bdg\n"
    cmd=cmd+"bedGraphToBigWig "+curr_sample+".sorted.bdg "+len_fn+" "+curr_sample+".bw\n"

    """
    Create UCSC track file
    """

    outp=open(track_fn, 'a')
    outp.write("track type=bigWig name="+"\""+curr_sample+"\" color="+curr_color+" gridDefault=on maxHeightPixels=50 visibility=full autoScale=off viewLimits=5:100 description=\""+curr_sample+"\" bigDataUrl="+bigdata_path+curr_sample+".bw\n")
    outp.close()

    return cmd

def lsf_file(job_name, cmd, memory=36000, thread=12):
    """
    Creates .lsf files
    """

    outp = open(job_name+".lsf",'w')
    outp.write("#!/bin/bash\n")
    outp.write("#BSUB -L /bin/bash\n")
    outp.write("#BSUB -J "+job_name+"\n")
    outp.write("#BSUB -q normal\n")
    outp.write("#BSUB -o "+job_name+"_%J.out\n")
    outp.write("#BSUB -e "+job_name+"_%J.screen\n")
    outp.write("#BSUB -M "+str(memory)+"\n")
    outp.write("#BSUB -n "+str(thread)+"\n")
    outp.write(cmd)
    outp.write("\n")
    outp.close()


def main(sample_info_file, project_name, ref_genome, library_type, index_type, path_start, template_dir, bam2bw):
    """
    Read in phenotype sample info provided by the users and perform the following steps:
    1) Perform adapter trimming - if Illumina indexes are not available in phenotype files
    2) run fastqc if they are not available
    3) Get unique reads 
    4) Run bwa to align reads to reference genome
    5) Obtain various QC metrics on aligned files
    """

    import chipseq_userdefine_variables as userdef # read in user-defined variable python script

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


    # Set up reference genome files
    fa, bwa_index_prefix = get_genome_ref_files(ref_genome)
    print "ref_genome = "+ref_genome
    print "bwa_index_prefix = "+ bwa_index_prefix

    # Obtain sample information from phenotype file
    print "sample_info_file = " + sample_info_file
    check_exist(sample_info_file)
    info_dict = get_sample_info(sample_info_file)
    sample_names = info_dict["Sample"]

    # Check if R1 and/or R2 column for .fastq file path exists
    print "library_type = "+library_type
    if library_type not in ['SE','PE']:
        print "Please specify 'SE' or 'PE' to library_type. Please check!"

    if "R1" not in info_dict:
        print "R1 column for Read 1 .fastq paths does not exist. Please check!"
        sys.exit()
    if library_type in ["PE"]:
        if "R2" not in info_dict:
            print "Paired-end library is specified. R2 column for Read 2 .fastq paths does not exist. Please check!"
            sys.exit()

    # Check if Index column exists and create index list
    # If index does not exist, assign NA
    if "Index" not in info_dict:
        print "No Index column in phenotype file. All samples skip adapter trimming."
        indexes = ['NA']*len(sample_names)
    else:
        indexes=info_dict['Index']

        if all(map(lambda x: 'NA' in x, indexes)): # if all indexes are NA
            print "All specified indexes are NA. Skip adapter trimming."

        # Obtain dictionary with adapter sequences of corresponding index type
        index_dict=index_check(indexes, index_type, template_dir)

    # If perform bam to bigwig file conversion, create url and colors
    if bam2bw:
        # overwrite ucsc track file if it already exists
        track_fn=path_start+project_name+"_ucsc_track.txt"
        if os.path.exists(track_fn):
            print "Warning: ucsc track file already exists: "+track_fn+". Overwrite it."
            outp=open(track_fn,'w')
            outp.write("")
            outp.close()

        # read in url
        bigdata_url=userdef.bigdata_url
        bigdata_path=bigdata_url+"/"+project_name+"/"
        print "Convert .bam files to .bw files. Use user-provided URL: "+bigdata_url
        print "Need to copy all the generated .bw files under this path: "+bigdata_path

        # obtain genome length file
        if ref_genome == "hg38":
    	    len_fn = userdef.hg38_len

        elif ref_genome == "hg19":
    	    len_fn = userdef.hg19_len

        check_exist(len_fn)

        # create color list for ucsc track display based on treatment condition
        rgb_colors=["27,158,119", "217,95,2", '117,112,179', '231,41,138', '102,166,30', '230,171,2', '166,118,29', '102,102,102'] # Dark2 color set "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666" convert to RGB (https://www.rapidtables.com/web/color/RGB_Color.html)
        if "Treatment" not in info_dict:
            colors=track_colors[0]*len(sample_names) # use one color for all samples
        else:
            conds=info_dict["Treatment"]
            cond_uniq=list(set(conds))
            # Assign color to each Status condition
            RGBS={}
            for i in range(len(cond_uniq)):
                RGBS[cond_uniq[i]]=rgb_colors[i]
            colors=map(lambda x: RGBS[x], conds)

    ####
    # Run by each sample
    ####

    for i in range(len(sample_names)):
        curr_sample=sample_names[i]
        curr_index=indexes[i]

        # Obtain read1 and read2 .fastq file paths
        R1=info_dict["R1"][i]
        check_exist(R1)
        if library_type in ["PE"]:
            R2=info_dict["R2"][i]
            check_exist(R2)
        elif library_type in ["SE"]:
            R2=""

        # Create output directory
        out_dir = path_start+curr_sample+"/"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Create cmd varialble to save linux commands output in .lsf files
        cmd = "cd "+out_dir+"\n"

        ###
        # Trim adaptor and fastqc
        ###
        if curr_index=="NA": # no index specified. skip trim
            R1_trim = R1
            R2_trim = R2
        else:
            # create adapter .fa file
            make_adapter_fa(curr_sample, out_dir, curr_index, index_dict)
            # create trimmed filename
            R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
            if R2=="":
                R2_trim = ""
            else:
                R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"
            
        # trim_and_fastqcL 1. trim, 2. fastqc, 3 total unique counts from .fastq file
        trim_and_fastqc_cmd=trim_and_fastqc(curr_sample, curr_index, out_dir, R1, R2, R1_trim, R2_trim)
        cmd = cmd + trim_and_fastqc_cmd

        ###
        # Align and obtain QC metrics
        ###

        # bwa alignment and bam stats
        bwa_cmd = bwa_and_bamstat(curr_sample, out_dir, library_type, R1_trim, R2_trim, bwa_index_prefix)
        cmd = cmd + bwa_cmd

        ###
        # Convert bam to bw
        ###

        if bam2bw:
            curr_color=colors[i]
            bam2bw_cmd=bam2bw_and_track(curr_sample, curr_color, out_dir, template_dir, track_fn, bigdata_path, len_fn)
            cmd=cmd+bam2bw_cmd

        ###
        # Create .lsf files
        ###
        lsf_file(curr_sample+"_align", cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for ChIP-seq samples associated with a project.")
    parser.add_argument("--project_name", type=str, help="Prefix name of for output directory and files")
    parser.add_argument("--samples_in", help="A tab-delimited txt file containing sample information with full path. See example file: sample_info_file.txt")
    parser.add_argument("--ref_genome", default="hg38", type=str, help="Specify reference genome (options: hg38, hg19, mm38, mm10, rn6, susScr3)")
    parser.add_argument("--library_type", default="PE", type=str, help="Specify library type (options: PE (paired-end), SE (single-end))")
    parser.add_argument("--index_type", type=str, help="If Index column is in phenotype file, specify index type for adapter trim."
        "(options: truseq_single_index, illumin_ud_sys1, illumin_ud_sys2 or user specified in the user-defined adapter reference file.)")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
    parser.add_argument("--template_dir", default="./", type=str, help="directory to put provided or user defined reference index files")
    parser.add_argument("--bam2bw", action='store_true', help="If specified, generate bigwig files (.bw) and create ucsc track file.")
    args = parser.parse_args()

    if args.project_name is None or args.samples_in is None:
        parser.print_help()
        sys.exit()

    main(args.samples_in, args.project_name, args.ref_genome, args.library_type, args.index_type, args.path_start, args.template_dir, args.bam2bw)

