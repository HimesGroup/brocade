#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import re

def get_sample_info(fin):
    """
    Read in information from phenotype file
    Create a directory to store each column
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

###
# Read in and process insert size files
###

def read_InsertSizeMetrics(fin):
	"""
	Read in output file created by Picardtools CollectInsertSizeMetrics function 
	Reformat metrics summary and histogram data to be put into table/plot in Rmd report
	"""
	f = open(fin,'r')
	c = f.read().split('\n\n')[1:]
	metrics = c[0].split('\n')[1:]
	hist = c[1].split('\n')[1:]
	metrics_out, hist_out = [], []
	for x in metrics:
		metrics_out.append( x.split('\t') )
	for x in hist:
		hist_out.append( x.split('\t') )
	return metrics_out, hist_out
	

def make_insertsize_matrix(path_out, project_name, sample_names, sample_paths, library_type):

	if library_type in ["PE"]:
		outp7 = open(path_out+project_name+"_insertmetrics_summary.txt", "w")
		name7 = ["Type"]
		for i in range(len(sample_names)):
			curr_name = sample_names[i]
			curr_path = sample_paths[i]
			name7.append(curr_name)
			#Make individual insert size files for each sample
			outp8 = open(path_out+project_name+"_"+curr_name+"_insertmetrics_hist.txt", "w")
			name8 = ["Insert_Size", curr_name]
			outp8.write("\t".join(name8)+"\n")
			insert_summary, insert_hist = read_InsertSizeMetrics(curr_path+"bwa_out/"+curr_name+".InsertSizeMetrics")
			if i == 0:
				g = zip(insert_summary[0], insert_summary[1])
				h = insert_hist
			else:
				for j in range(len(g)):
					g[j] = list(g[j])+[insert_summary[1][j]]
				h = insert_hist
			outp8.write("\n".join(map("\t".join, h[1:])))
			outp8.close()
		outp7.write("\t".join(name7)+"\n")
		outp7.write("\n".join(map("\t".join, g) ))
		outp7.close()
		print "Created file "+path_out+project_name+"_insertmetrics_summary.txt and sample-specific *_insertmetrics_hist.txt files"

###
# Read in and process samtools idxstat files
###


def read_samtools_stats(fin, ref_genome):
    """
    Read in output file created by Samtools stats function (of type curr_sample.stats)
    Reformat and output:
    ref_genome summary data to be put into table and plot in Rmd report
    """
    f = open(fin,'r')
    c = f.read().split('\n')
    if '' in c:
        c.remove('')
    rna_out = []
    rrna, other = 0, 0
    for x in c:
	if ref_genome == "hg19":
	    if ("chrM" in x) or ("chrUn_gl000220" in x):
		rrna += int(x.split('\t')[2])
	    elif len(x.split('\t')[0].split('_')) > 1:
		other += int(x.split('\t')[2])
	    elif "*" not in x:
		curr_line = x.split('\t')[:-1]
		curr_line[0] = curr_line[0].strip("chr")
                rna_out.append( curr_line )
	
	if ref_genome == "hg38":
	    if ("chrM" in x) or ("chrUn_GL000220v1" in x):
		rrna += int(x.split('\t')[2])
	    elif len(x.split('\t')[0].split('_')) > 1:
		other += int(x.split('\t')[2])
	    elif ("*" not in x) and ("EBV" not in x):
		curr_line = x.split('\t')[:-1]
		curr_line[0] = curr_line[0].strip("chr")
		rna_out.append( curr_line )

	if ref_genome == "mm10" or ref_genome == "mm38":
	    if len(x.split('\t')[0].split('.')) > 1:
		other += int(x.split('\t')[2])
	    elif len(x.split('\t')[0].split('_')) > 1:
		other += int(x.split('\t')[2])
	    elif ".1" and "_" and "*" not in x:
		curr_line = x.split('\t')[:-1]
		curr_line[0] = curr_line[0]
		rna_out.append( curr_line )
	    rrna = "NA"

        if ref_genome == "rn6":
	    if len(x.split('\t')[0].split('.')) > 1:
		other += int(x.split('\t')[2])
	    elif len(x.split('\t')[0].split('_')) > 1:
		other += int(x.split('\t')[2])
	    elif ".1" and "_" and "*" not in x:
		curr_line = x.split('\t')[:-1]
		curr_line[0] = curr_line[0]
		rna_out.append( curr_line )
            rrna = "NA"

	if ref_genome == "susScr3":
	    if ("chrM" in x) or ("chrUn_gl000220" in x):
		rrna += int(x.split('\t')[2])
	    elif len(x.split('\t')[0].split('_')) > 1:
		other += int(x.split('\t')[2])
	    elif "*" not in x:
		curr_line = x.split('\t')[:-1]
		curr_line[0] = curr_line[0].strip("chr")
		rna_out.append( curr_line )

	if ref_genome == "Zv9":
	    if len(x.split('\t')[0].split('_')) > 1:
		other += int(x.split('\t')[2])
	    elif "*" not in x:
		curr_line = x.split('\t')[:-1]
		curr_line[0] = curr_line[0]
		rna_out.append( curr_line )
	    rrna = "NA"
    rna_out.append(['Other','',str(other)])
    rna_out.append(['rRNA','',str(rrna)])

    return rna_out



def make_samidxstat_matrix(path_out, project_name, sample_names, sample_paths, ref_genome):
    """
    Read in individual output files created by curr_path+"/"+bwa_out+"/"+curr_name+".stats" according to sample_info_file
    Return as a single matrix for a whole batch
    """

    # path_out: directory for output files
    # project_name: file prefix
    # sample_names: a list stores all sample names
    # sample_paths: a list stores all directory paths for each sample
    # ref_genome: reference genome

    name1 = ["Chromosome", "Length"] # create file header
    c=[] # a list to store samtools idx stats for each sample by chromosome

    for i in range(len(sample_names)):
        curr_name = sample_names[i]
	curr_path = sample_paths[i]+"bwa_out/"+curr_name+".stats"
	if not os.path.exists(curr_path):
            print "Missing samtools idxstat output file ", curr_path
	    break
 
        rna_out = read_samtools_stats(curr_path, ref_genome)
        #print "Read in samtools idxstat output for sample "+curr_name+": Done"
        name1.append(curr_name)

	if i == 0:
            c = rna_out # for the first sample
	else:
	    for j in range(len(c)):
                c[j] = c[j]+[rna_out[j][2]]

        #print "Append sample "+curr_name+" to matrix: Done"

    # output samtools idx stat file
    outp1 = open(path_out+project_name+"_counts.txt", "w")
    outp1.write("\t".join(name1)+"\n")
    outp1.write("\n".join(map("\t".join, c)))
    outp1.write("\n")
    outp1.close()
    print "Created samtools idx stat matrix file "+path_out+project_name+"_counts.txt"


###
# Read in and process bamstats files
###

def read_bamtools_stats(fin):
	"""
	Read in bamstats output file curr_sample.bamstats
	Reformats and outputs portions of interest for Rmd report
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[5:]
	for x in range(c.count('')):
		c.remove('')
	bs = map(lambda x: x.split(':'), c)
	names = map(lambda x: x[0], bs)
	counts = map(lambda x: x[1].split('\t')[0].strip(' '), bs)
	bamstats = zip(names, counts)
	return bamstats

def make_bamstats_matrix(path_out, project_name, sample_names, sample_paths):
	outp6 = open(path_out+project_name+"_bamstats_counts.txt", "w")
	name6 = ["Type"]
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		name6.append(curr_name)
		bam_stats = read_bamtools_stats(curr_path+"bwa_out/"+curr_name+".bamstats")
		if i == 0:
			f = bam_stats
		else:
			for j in range(len(bam_stats)):
				f[j] = list(f[j])+[bam_stats[j][1]]
	outp6.write("\t".join(name6)+"\n")
	outp6.write("\n".join(map("\t".join, f)))
        outp6.write("\n")
	outp6.close()
	print "Created file "+path_out+project_name+"_bamstats_counts.txt"


###
# Read in and process _ReadCount from raw .fastq files
###

def get_unique_reads(fin, library_type):
	"""
	Read in file that has output from awk on number of reads, unique reads, and %unique reads from fastq files
	Output R1 number of reads, R1 unique reads, R1 %unique reads, R2 number of reads, R2 unique reads, and R2 %unique reads for Rmd report
	"""
	f = open(fin,'r')
	c = f.read().split('\n')
	if '' in c:
		c.remove('')
	#First row is R1. Second row is R2. 
	read_numbers = map(lambda x: x.split(' '), c)
	if library_type in ["PE"]:
		read_numbers = read_numbers[0]+read_numbers[1]
	else:
		read_numbers = read_numbers[0]
	return read_numbers

def make_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type):
    outp9 = open(path_out+project_name+"_unique_counts.txt", "w")
    if library_type in ["PE"]:
        name9 = ["Sample", "R1_Raw_Reads", "R1_Unique_Reads", "R1_Percent_Unique", "R2_Raw_Reads", "R2_Unique_Reads", "R2_Percent_Unique"]
    else:
	name9 = ["Sample", "Raw_Reads", "Unique_Reads", "Percent_Unique"]
    for i in range(len(sample_names)):
	curr_name = sample_names[i]
	curr_path = sample_paths[i]
	unique_reads = get_unique_reads(curr_path+curr_name+"_ReadCount", library_type)
	if i == 0:
            unique_counts = [[curr_name] + unique_reads]
	else:
	    unique_counts.append( [curr_name] + unique_reads)

    outp9.write("\t".join(name9)+"\n")
    outp9.write("\n".join(map("\t".join, unique_counts)))
    outp9.write("\n")
    outp9.close()
    print "Created file "+path_out+project_name+"_unique_counts.txt"


###
# Read in and process fastqc files
###

def read_fastq_data(fin):
	"""
	Read fastqc_data.txt file from a FastQC report and extract the percentage of duplicates and histogram information
	"""
	f = open(fin,'r')
	c = f.read().split('#Total Deduplicated Percentage')[1]
	c2 = c.split('>>END_MODULE')[0]
	c3 = c2.split('\n')
	if '' in c3:
		c3.remove('')
	#Each entry contains Duplication Level, Percentage of deduplicated, Percentage of total
	duplicate_data =  [['Total Deduplicated Percentage', c3[0].strip('\t')]]+map(lambda x: x.split('\t'), c3[2:])
	return duplicate_data

def make_duplicate_matrix(path_out, project_name, sample_names, sample_paths, library_type):

    #Duplicate read info from fastq files
    outp10 = open(path_out+project_name+"_duplicates.txt", "w")
    name10 = ["Read_Number"]
    for i in range(len(sample_names)):
        curr_name = sample_names[i]
        curr_path = sample_paths[i]
	#Move individual FastQC reports to report directory
	subprocess.call("cp "+curr_path+curr_name+"_R1_Trimmed_fastqc.zip "+path_out, shell=True)
	subprocess.call("unzip -o -q -d "+path_out+" "+path_out+curr_name+"_R1_Trimmed_fastqc.zip", shell=True)
	subprocess.call("rm  "+path_out+curr_name+"_R1_Trimmed_fastqc.zip", shell=True)
	name10.append(curr_name+"_R1")	 
	if os.path.isfile(path_out+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"):
	    fastqc1 = path_out+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"
	elif os.path.isfile(path_out+"/"+curr_name+"_R1_fastqc/fastqc_data.txt"):
	    fastqc1 = path_out+"/"+curr_name+"_R1_fastqc/fastqc_data.txt"
	elif os.path.isfile(path_out+"/"+curr_name+"_fastqc/fastqc_data.txt"):
	    fastqc1 = path_out+"/"+curr_name+"_fastqc/fastqc_data.txt"
	else:
	    print "Missing FastQC report", curr_name
	    break

	curr_dup_R1 = read_fastq_data(fastqc1)

	if library_type in ["PE"]:
	    subprocess.call("cp "+curr_path+curr_name+"_R2_Trimmed_fastqc.zip "+path_out, shell=True)
	    subprocess.call("unzip -o -q -d "+path_out+" "+path_out+curr_name+"_R2_Trimmed_fastqc.zip", shell=True)
	    subprocess.call("rm  "+path_out+curr_name+"_R2_Trimmed_fastqc.zip", shell=True)
	    name10.append(curr_name+"_R2")
	    if os.path.isfile(path_out+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"):
                fastqc2 = path_out+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"
	    elif os.path.isfile(path_out+"/"+curr_name+"_R2_fastqc/fastqc_data.txt"):
		fastqc2 = path_out+"/"+curr_name+"_R2_fastqc/fastqc_data.txt"
	    else:
		print "Missing FastQC report for R2 ", curr_name
	    curr_dup_R2 = read_fastq_data(fastqc2)

	if i == 0:
	    duplicates = curr_dup_R1
	    for j in range(1, len(duplicates)):
                duplicates[j] = [duplicates[j][0], duplicates[j][2]]
	    if library_type in ["PE"]:
                duplicates[0].append(curr_dup_R2[0][1])
		for j in range(1, len(duplicates)):
                    duplicates[j].append(curr_dup_R2[j][2])
	else:
	    duplicates[0].append(curr_dup_R1[0][1])
	    for j in range(1, len(duplicates)):
	        duplicates[j].append(curr_dup_R1[j][2])
	    if library_type in ["PE"]:
	        duplicates[0].append(curr_dup_R2[0][1])
		for j in range(1, len(duplicates)):
                    duplicates[j].append(curr_dup_R2[j][2])

    outp10.write("\t".join(name10)+"\n")
    outp10.write("\n".join(map("\t".join, duplicates)))
    outp10.write("\n")
    outp10.close()
    print "Created file "+path_out+project_name+"_duplicates.txt"


def make_project_data_files(project_name, sample_names, sample_paths, path_out, ref_genome, library_type):
	"""
	Creates several text files to be loaded by R for Rmd report based on modified program outputs read in with preceding scripts
	These text files are matrices containing information for all samples in a project/batch
	Currently, cycles through all samples multiple times, once to create each individual file type
	Currently, there is no way to handle missing files. If an error is encountered the process will stop at that point and not complete
	Currently, assumes default naming convention of all programs used in chipseq_align_and_qc.py
	"""

	#Unique read counts obtained by comprehensive count of fastq files -- outp9
        make_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type)
	
	#Duplicate read info from fastq files -- outp10
        make_duplicate_matrix(path_out, project_name, sample_names, sample_paths, library_type)

	#samtools stats counts of reads per chromosome
        make_samidxstat_matrix(path_out, project_name, sample_names, sample_paths, ref_genome)

	#bamstats output metrics on types of reads, including junction spanning reads
        make_bamstats_matrix(path_out, project_name, sample_names, sample_paths)

	#insertsizemetrics output on insert size statistics -- outp6
        make_insertsize_matrix(path_out, project_name, sample_names, sample_paths, library_type)

def make_rmd_css(new_dir):
    """
    create custom.css for rmarkdown
    """

    css_outp = open(new_dir+"custom.css", "w")
    css_outp.write("blockquote {\n")
    css_outp.write("    padding: 10px 20px;\n")
    css_outp.write("    margin: 0 0 20px;\n")
    css_outp.write("    font-size: 14px;\n")
    css_outp.write("    background-color: #eee;\n")
    css_outp.write("    border-left: 5px solid #eee;\n")
    css_outp.write("}\n\n")
    css_outp.write(".main-container {\n")
    css_outp.write("    max-width: 2000px !important;\n") # "!important" overrides other rules
    css_outp.write("}\n")
    css_outp.close()


def make_rmd_title(new_dir, project_name, ref_genome):
    """
    Rmarkdown creation: Create title and major description
    """

    import chipseq_userdefine_variables as userdef # read in user-defined variable python script
    # import software version
    trimmomatic_version=userdef.trimmomatic_version
    fastqc_version=userdef.fastqc_version
    bwa_version=userdef.bwa_version
    samtools_version=userdef.samtools_version
    bamtools_version=userdef.bamtools_version
    picard_version=userdef.picard_version
    # import author information
    author=userdef.author

    rmd="---\ntitle: 'ChIP-Seq Report of Sample QC and Alignment Summary Statistics for "+project_name+"'\n"
    rmd=rmd+"author: "+author+"\n"
    rmd=rmd+"date: \"`r format(Sys.time(), '%d %B, %Y')`\"\n"
    rmd=rmd+"output: \n"
    rmd=rmd+"  html_document:\n"
    rmd=rmd+"    css: custom.css\n"
    rmd=rmd+"    toc: true\n"
    rmd=rmd+"    toc_float: true\n---\n\n"

    # Variable used
    rmd=rmd+"**Project:** "+project_name+"\n\n"
    if ref_genome == "hg19":
        rmd=rmd+"**Genome:** For human, the hg19 assembly was used. We estimate the number of rRNA reads as those mapped to chrM plus chrUn_gl000220, corresponding to 12S, 16S and 5.8S rRNA. The 'Other' category contains all other chr*_random and chrUn_* available. If using the 2014 updated version of the hg19 files, these categories are no longer present.\n"
    elif ref_genome == "hg38":
        rmd=rmd+"**Genome:** For human, the hg38 assembly was used. We estimate the number of rRNA reads as those mapped to chrM plus chrUn_GL000220v1, corresponding to 12S, 16S and 5.8S rRNA. The 'Other' category contains all other chr*_random and chrUn_* available.\n"
    elif ref_genome == "mm38":
	rmd=rmd+"**Genome:** For mouse, the ENSEMBL GRCm38 assembly available in iGenomes was used.\n"
    elif ref_genome == "mm10":
	rmd=rmd+"**Genome:** For mouse, the UCSC mm10 assembly available in iGenomes was used.\n"
    elif ref_genome == "rn6":
	rmd=rmd+"**Genome:** For rat, the rn6 assembly was used.\n"
    elif ref_genome == "susScr3":
	rmd=rmd+"**Genome:** For pig, the susScr3 assembly was used.\n"
    elif ref_genome == "Zv9":
	rmd=rmd+"**Genome:** For zebrafish, the Zv9 assembly comprises a sequence length of 1.4 Gb in 26 chromosomes (labels 1-25 and MT) and 1,107 scaffolds (merged into label 'Other').\n"
    rmd=rmd+"\n\n"

    # Bioinformatics tools
    rmd=rmd+"**Informatics tools used:**\n\n"
    rmd=rmd+"* Trimmomatic ("+trimmomatic_version+")\n"
    rmd=rmd+"* FastQC ("+fastqc_version+")\n"
    rmd=rmd+"* BWA ("+bwa_version+")\n"
    rmd=rmd+"* samtools ("+samtools_version+")\n"
    rmd=rmd+"* bamtools ("+bamtools_version+")\n"
    rmd=rmd+"* Picard Tools ("+picard_version+")\n"
    rmd=rmd+"\n\n"

    return rmd

def make_rmd_var(project_name, path_start, ref_genome, library_type, sample_names, sample_info_file, template_dir):
    """
    Rmarkdown creation: Define variables and files
    """

    rmd="```{r vars, echo=F}\n"
    rmd=rmd+"project_name=\""+project_name+"\"\n"
    rmd=rmd+"path_start=\""+path_start+"\"\n"
    rmd=rmd+"sample_names <- c("+str(sample_names)[1:-1]+")\n"
    rmd=rmd+"genome=\""+ref_genome+"\"\n"
    rmd=rmd+"library_type=\""+library_type+"\"\n"
    rmd=rmd+"sample_info_file='"+sample_info_file+"'\n"
    rmd=rmd+"```\n\n"

    return rmd



def make_rmd_featurestat(project_name):
    """
    Rmarkdown creation: Obtain proportion of no feature counts from htseq-count results
    Separate here because of the htseq-count-specific feature
    """
    rmd="\n"
    rmd=rmd+"## HTSeq-count: No feature counts statistics\n\n"
    rmd=rmd+"No feature count (per million reads) statistics from htseq-count quantification results\n\n"
    rmd=rmd+"```{r nofeature, eval=T, echo=F, message=F, warning=F, results='asis'}\n"
    rmd=rmd+"nofeature.data <- read.table('"+project_name+"_htseq_nofeature.txt', sep='\\t', header=T, as.is=T)\n"
    rmd=rmd+"DT::datatable(nofeature.data)\n"
    rmd=rmd+"```\n\n"
    return rmd

def lsf_file(job_name, cmd, memory=36000, thread=1):
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

def make_rmd_html(sample_info_file, project_name, path_start, new_dir, sample_names, ref_genome, library_type, template_dir):
    """
    Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
    Runs rmarkdown to create html document
    Also makes a custom css format file - improvement on rmarkdown defaults
    """
    # create custom.css for rmarkdown
    make_rmd_css(new_dir)

    # create contents for rmarkdown RMD file
    rmd = ""

    # create title and description
    rmd_title =make_rmd_title(new_dir, project_name, ref_genome)
    rmd = rmd + rmd_title

    # create varialbes
    rmd_var = make_rmd_var(project_name, path_start, ref_genome, library_type, sample_names, sample_info_file, template_dir)
    rmd = rmd + rmd_var

    # read in contents in template Rmdfile
    rmd_in = open(template_dir+"chipseq_align_and_qc_report_Rmd_template.txt", "r")
    rmd_template = rmd_in.read()
    rmd = rmd + rmd_template
    rmd = rmd + "\n\n"

    # output no feature count stat table if htseq-count is used
    #if aligner=="star":
    #   rmd_featurestat = make_rmd_featurestat(project_name)
    #   rmd = rmd + rmd_featurestat

    # write in project_name+"_QC_ChIPSeqReport.Rmd" file
    outp = open(new_dir+project_name+"_QC_ChIPSeqReport.Rmd", "w")
    outp.write(rmd)

    outp.write("\n\n")
    outp.close()
    print "Created file "+new_dir+project_name+"_QC_ChIPSeqReport.Rmd"

    # create .lsf file. Run it on HPC becuase rlog for all sample counts in pca step is computationally demanding.
    cmd="cd "+new_dir+"; echo \"library(rmarkdown); rmarkdown::render('"+project_name+"_QC_ChIPSeqReport.Rmd')\" | R --no-save --no-restore\n"
    lsf_file(project_name+"_qc", cmd)
    print "Created LSF script "+project_name+"_qc.lsf in current directory"

def main(project_name, sample_info_file, path_start, ref_genome, library_type, template_dir):
    """
    Creates html report describing summary and QC statistics for a set of aligned ChIP-Seq samples associated with a project
    Report is based on multiple output files created by chipseq_align_and_qc.py
    Such files are first reformatted into matrices that are easily loaded into R
    Input:
        project_name: name for report.
	sample_info_file: tab delimited txt file with sample information as described in chipseq_align_and_qc.py
	rmd_template: txt file that contains most of the contents that will populate the Rmd file for the report
    Current ref_genome choices: hg19, mm38
    Current library_type choices: PE, SE, DGE, SPE
    """
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
	path_start = path_start+"/"
    new_dir = path_start+project_name+"_Alignment_QC_Report/"
    if not os.path.exists(new_dir):
	os.makedirs(new_dir)
    if template_dir == "./":
        template_dir = os.getcwd()
    if template_dir[-1] != "/":
        template_dir = template_dir+"/"

    # check if QC template txt file exists
    if not os.path.exists(template_dir+"chipseq_align_and_qc_report_Rmd_template.txt"):
        print "Cannot find "+template_dir+"chipseq_align_and_qc_report_Rmd_template.txt"
	sys.exit()

    #Get list of dictionary of sample information. Keys: [sample_id, ercc_mix, index]
    info_dict = get_sample_info(sample_info_file)
    sample_names = info_dict["Sample"]
    sample_paths = map(lambda x: path_start+x+"/", sample_names)

    # create stats files used for QC report
    make_project_data_files(project_name, sample_names, sample_paths, new_dir, ref_genome, library_type)

    #Create the report
    make_rmd_html(sample_info_file, project_name, path_start, new_dir, sample_names, ref_genome, library_type, template_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create HTML report of QC and alignment summary statistics for ChIP-seq samples associated with a project.")
    parser.add_argument("--project_name", type=str, help="Prefix name of for output directory and files")
    parser.add_argument("--samples_in", help="A tab-delimited txt file containing sample information with full path. See example file: sample_info_file.txt")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
    parser.add_argument("--ref_genome", default="hg38", type=str, help="Specify reference genome (options: hg38, hg19, mm38, mm10, rn6, susScr3)")
    parser.add_argument("--library_type", default="PE", type=str, help="Specify library type (options: PE (paired-end), SE (single-end))")
    parser.add_argument("--template_dir", default="./", type=str, help="directory to put template RMD file seq_align_and_qc_report_Rmd_template.txt for QC report")
    args = parser.parse_args()

    if args.project_name is None or args.samples_in is None:
        parser.print_help()
        sys.exit()

    main(args.project_name, args.samples_in, args.path_start, args.ref_genome, args.library_type, args.template_dir)
