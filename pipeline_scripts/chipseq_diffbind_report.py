#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import re
import fnmatch

# import user-defined parameters
import chipseq_userdefine_variables as userdef # read in user-defined variable python script
# import software version
bwa_version=userdef.bwa_version
macs2_version=userdef.macs2_version
diffbind_version=userdef.diffbind_version
chipseeker_version=userdef.chipseeker_version
deseq2_version=userdef.deseq2_version
# import author information
author=userdef.author
# import favorite genes
fav_gene=userdef.fav_gene
# import HPC parameters
memory=userdef.memory
queue=userdef.queue

def make_rmd_css(path_start):
    """
    create custom.css for rmarkdown
    """

    css_outp = open(path_start+"custom.css", "w")
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


def make_diffbind_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file):
    """
    Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
    """

    ###
    # Create global design variable based on comparison file
    ###

    #load text file containing all comparisons of interest
    comps_file = open(comp_file)
    comps = comps_file.readlines()[1:] # exclude header line

    # Create out directory
    out_dir = path_start+project_name+"_diffbind_out/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ###
    # Create RMD file for diffbind
    ###

    # create custom.css for rmarkdown
    make_rmd_css(out_dir)

    outp = open(out_dir+'/'+project_name+"_DiffBind_Report.Rmd", "w")

    # title
    outp.write("---\ntitle: 'Differential Binding Site Results for "+project_name+"'\n")
    outp.write("author: "+author+"\n")
    outp.write("date: \"`r format(Sys.time(), '%d %B, %Y')`\"\n")
    outp.write("output: \n")
    outp.write("  html_document:\n")
    outp.write("    css: custom.css\n")
    outp.write("    toc: true\n")
    outp.write("    toc_float: true\n---\n\n")

    # description
    outp.write("Reads were aligned to the "+ref_genome+" assembly using BWA ("+bwa_version+").  The following alignment QC report was produced:<br>\n\n")
    outp.write("> "+project_name+"_DiffBind_Report.html<br>\n\n")
    outp.write("MACS2 ("+macs2_version+") function callpeak was used to call peaks. Peaks within blacklisted regions were filtered out.\n\n")
    outp.write("DiffBind ("+diffbind_version+") was used for differential binding site analaysis, based on the the phenotype file provided. Normalized counts are saved in the following text file:<br>\n\n")
    outp.write("> "+project_name+"_Condition1_vs_Condition0_counts_normalized_by_diffbind.csv<br>\n\n")
    outp.write("Differential binding site analysis was done for all comparisons provided in the comparisons file, using DESeq2 by default.<br>\n\n")
    if ref_genome=="hg38":
        txdb="TxDb.Hsapiens.UCSC.hg38.knownGene"
        annoDb="org.Hs.eg.db"
    elif ref_genome=="hg19":
        txdb="TxDb.Hsapiens.UCSC.hg19.knownGene"
        annoDb="org.Hs.eg.db"
    outp.write("Peaks were annotated using corresponding R packages "+ annoDb +" and "+ txdb +" by ChIPseeker ("+chipseeker_version+").<br>\n\n")

    # Load library and set variables
    outp.write("\n\n```{r lib, echo=F, message=F, warning=F}\n")
    outp.write("annoDb_name='"+annoDb+"'\n")
    outp.write("txdb_name='"+txdb+"'\n")
    outp.write("library(annoDb_name,character.only = TRUE)\nlibrary(txdb_name,character.only = TRUE)\n")
    outp.write("txdb=get(txdb_name)\n")
    outp.write("library(DiffBind)\nlibrary(ChIPseeker)\nlibrary(tidyr)\nlibrary(DT)\nlibrary(devtools)\nlibrary(ggplot2)\nlibrary(gplots)\nlibrary(RColorBrewer)\nlibrary(viridis)\nlibrary(pander)\noptions(width = 1000)\n```\n")
    outp.write("\n")
    outp.write("```{r vars, eval=T, echo=F}\n")
    outp.write("project_name=\""+project_name+"\"\n")
    outp.write("path_start='"+path_start+"'\n")
    outp.write("out_dir='"+out_dir+"'\n")
    outp.write("coldata <- read.table('"+sample_info_file+"', sep='\\t', header=TRUE)\n")
    outp.write("coldata <- subset(coldata, QC_Pass==1)\n")
    outp.write("tss_span=3000 # region range of transcription starting site\n")
    outp.write("flank_span=10000 # # region with peak +/-span\n")
    outp.write("```\n\n")

    # create DiffBind bam read files
    outp.write("Create bam read columns\n")
    outp.write("```{r bamReads, eval=T, echo=F}\n")
    outp.write("coldata$bamReads=paste0(path_start,coldata$Sample,'/bwa_out/',coldata$Sample,'.bam')\n")
    outp.write("bamreads=coldata$bamReads[!is.na(coldata$Input)]\n")
    outp.write("nobamReads=bamreads[!sapply(bamreads,file.exists)]\n")
    outp.write("if (length(nobamReads)>1) {stop('Bam read file(s) do not exists: ',paste(nobamReads,collapse=', '))}\n")
    outp.write("```\n\n")

    # create DiffBind bam control files
    outp.write("Create bam control columns\n")
    outp.write("```{r bamControl, eval=T, echo=F}\n")
    outp.write("coldata$bamControl=paste0(path_start,coldata$Input,'/bwa_out/',coldata$Input,'.bam')\n")
    outp.write("bamcontrols=coldata$bamControl[!is.na(coldata$Input)]\n")
    outp.write("nobamControl=bamcontrols[!sapply(bamcontrols,file.exists)]\n")
    outp.write("if (length(nobamControl)>1) {stop('Bam control file(s) do not exists: ',paste(nobamControl,collapse=', '))}\n")
    outp.write("```\n\n")

    # create DiffBind peak bed files
    outp.write("Create peak bed columns\n")
    outp.write("```{r bed, eval=T, echo=F}\n")
    outp.write("coldata$Peaks=paste0(path_start,coldata$Sample,'/macs2_out/',coldata$Sample,'.blackfilt.bed')\n")
    outp.write("bedpeaks=coldata$Peaks[!is.na(coldata$Input)]\n")
    outp.write("nobed=bedpeaks[!sapply(bamcontrols,file.exists)]\n")
    outp.write("if (length(nobed)>1) {stop('Peak bed file(s) do not exists: ',paste(nobed,collapse=', '))}\n")
    outp.write("```\n\n")



    #create and paste the portion of the report that is unique to each comparison
    for line in comps:
        line=line.rstrip()
        case=line.split('\t')[0]
        ctrl=line.split('\t')[1]

        outp.write("## "+case+" vs. "+ctrl+"\n")
        outp.write("\n")
        outp.write("```{r, eval=T, echo=F}\n")
        outp.write("cond1 <- '"+case+"'\n")
        outp.write("cond0 <- '"+ctrl+"'\n")
        #outp.write("res <- results(dds, contrast=c('Status','"+case+"','"+ctrl+"'))\n")
        outp.write("```\n\n")

	outp.writelines(rmd_template)
	outp.write("\n")

    outp.write("```{r session_info, eval=T, echo=F}\n")
    outp.write("pander(sessionInfo())\n")
    outp.write("```\n\n")
    outp.close()

    # create .lsf file for HPC use
    lsf_cmd="cd "+out_dir+"; echo \"library(rmarkdown); rmarkdown::render('"+project_name+"_DiffBind_Report.Rmd')\" | R --no-save --no-restore\n"
    lsf_file(project_name+"_diffbind", lsf_cmd)


def main(project_name, sample_info_file, path_start, comp_file, template_dir, ref_genome):
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
        path_start = path_start+"/"

    if template_dir == "./":
	template_dir = os.getcwd()
    if template_dir[-1] != "/":
        template_dir = template_dir+"/"

    # check if sample info exists
    if not os.path.exists(sample_info_file):
        print "Cannot find sample_info_file: "+sample_info_file
        sys.exit()

    # check if compare file exists
    if not os.path.exists(comp_file):
        print "Cannot find the comparison file: "+comp_file
        sys.exit()

    # check if diffbine template file exists
    if not os.path.exists(template_dir+"chipseq_diffbind_Rmd_template.txt"):
        print "Cannot find chipseq_diffbind_Rmd_template.txt"
	sys.exit()

    rmd_in = open(template_dir+"chipseq_diffbind_Rmd_template.txt", "r")
    rmd_template = rmd_in.readlines()
    make_diffbind_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create HTML report of differential expression results for ChIP-seq samples associated with a project.")
    parser.add_argument("--project_name", type=str, help="Prefix name of project for all output files.")
    parser.add_argument("--samples_in", help="A tab-delimited txt file containing sample information with full path. See example file: sample_info_file.txt, but add an additional QC_Pass column")
    parser.add_argument("--comp", help="A tab-delimited txt file containing sample comparisons to be made. One comparison per line, columns are Condition1, Condition0, Design. "
            "Design: specify paired or unpaired. For paired design, specify condition to correct for, matching the column name in the 'coldata' file - e.g. paired:Donor.")
    parser.add_argument("--ref_genome", default="hg38", type=str, help="Specify reference genome (options: hg38, hg19)")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path to where project directory created by chipseq_diff_report.py is located (default=./)")
    parser.add_argument("--template_dir", default="./", type=str, help="directory to put template RMD file chipseq_diffbind_Rmd_template.txt for QC report")

    args = parser.parse_args()

    if args.comp is None or args.project_name is None or args.samples_in is None:
        parser.print_help()
        sys.exit()

    main(args.project_name, args.samples_in, args.path_start, args.comp, args.template_dir, args.ref_genome)


