#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Jun  6 14:42:00 2018
Author: Avantika Diwadkar (@diwadkar)


This script can be used to generate a phenotype file for any given ChIP-Seq dataset.
This is followed by using the this information to acquire SRR numbers and download the 
fastq files from the server. 

"""

#Libraries 
import os
import subprocess
import urllib
import urllib.request



#Functions used
def make_dir(filepath):
    if not os.path.exists(filepath):
        os.makedirs(filepath)
        print(filepath.split("/")[-1] + " directory created in " + filepath)
    

def download_file(url, filepath):
    if not os.path.isfile(filepath):
        urllib.request.urlretrieve(url, filepath)
        print("Downloaded " + filepath.split("/")[-1] + "file in " + filepath)
        
        
def write_in_file(file, codestring):
    pfile = open(file,"w")
    pfile.write(codestring)
    pfile.close()
    

#Acquire matrix file for phenotype information

#Take input geo_id to set variables
geo_id = input("Please enter GEO ID of ChIP-Seq dataset: ")
matrix_file_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/" + geo_id[:-3] + "nnn/" + geo_id + "/matrix/"


#Set path for download
path = os.getcwd() + "/"
geo_dir = path+geo_id
data_dir = geo_dir + "/" + "data"
res_dir = geo_dir + "/" + "results"

    

#Download matrix file in data directory if it does not exist
file_path = data_dir+ "/" + geo_id + "_series_matrix.txt.gz"

#As urllib request is giving corrupted file, use wget
matrix_file = matrix_file_url + geo_id + "_series_matrix.txt.gz"
mat_cmd = "wget %s"%(matrix_file)


#Make chipseq analysis file: edit phenotype information

Rcall = '''
---
title: "SRA downlaod for brocade"
output: html_document
---

```{r, echo=FALSE}
knitr::opts_chunk$set(error = TRUE)
```
    
```{r setup, eval=F, echo=F, message=F, warning=F}
#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite ("SRAdb")
```

Loading all the essential libraries
```{r lib, eval=T, echo=T, message=F, warning=F}
library(GEOquery)
library(SRAdb)
library(plyr)
library(dplyr)
library(pander)
```

#### Getting the ExpressionSet data from the matrix file
Read in the matrix data file and assign data and result folders
```{r gse, eval=T, echo=T, message=F, warning=F}
datadir = "data_dir"
resdir = "res_dir"
geo_id = "demo"
geo_fn = paste0(geo_id,"_series_matrix.txt.gz")
gse <- getGEO(filename=paste0(datadir,geo_fn),GSEMatrix = TRUE)
```

#### Phenotype file preparation

A) Read in the phenotype from the ExpressionSet
```{r pheno.raw, eval=T, echo=T, message=F, warning=F}
pheno.raw <- pData(phenoData(gse))
pandoc.table(sapply(pheno.raw,levels),caption="All phenotype levels in dataset")

```

B) Mutate and save phenotype information in text file
```{r pheno, eval=T, echo=T, message=F, warning=F}
pheno <- pheno.raw %>% select(title, geo_accession, organism_ch1, characteristics_ch1, characteristics_ch1.1, data_processing.4)
pheno <- as.data.frame(lapply(pheno, function(y) gsub(".*:", "", y)))
colnames(pheno) <- c("Sample","ID","Organism","Tissue","Treatment","Genome")
pheno$Antibody  <- gsub(" .*","", pheno$Sample)
pheno$Subject  <- gsub("_.*","", pheno$Sample) #Need to edit input type
write.table(pheno,paste0(resdir,geo_id,"_Phenotype_withoutQC.txt"),col.names=T,row.names=F,sep="\t",quote=F)
pandoc.table(pheno,split.tables=Inf, caption="Selected phenotype of the dataset")
```


#### Download the RAW data files from SRA archive using SRAdb package. 
Get SRA files
```{r sra, eval=T, echo=T, message=F, warning=F}
sqlfile <<- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
files <- gsub(".*=","",pheno.raw$relation.1)
groups <- gsub(".*:","",pheno.raw$characteristics_ch1.1)
dir.create("fastq",showWarnings = FALSE)
getSRAfile(files, sra_con,
  destDir = file.path(getwd(),"fastq") , fileType = 'fastq',
  srcType= 'ftp', makeDirectory = FALSE,
  method = 'libcurl', ascpCMD = NULL )
```


#### Session information

```{r sessioninfo, eval=T, echo=F}
pander(sessionInfo())
```

'''

#Make RMD file to prepare phenotype file and download SRA data
pheno_file = geo_dir + "/" + geo_id + "_analysis.Rmd" 

         
#Change phenotype file geo_id variable in chip-seq RMD
#For OS X add '' -e after -i
cmd = "sed -i 's/demo/%s/g' %s"%(geo_id,pheno_file)

#Change paths in RMD file
dat_cmd = "sed -i 's/datadir/%s/g' %s"%(data_dir,pheno_file)
res_cmd = "sed -i 's/resdir/%s/g' %s"%(res_dir,pheno_file)


#Make RMD markdown file - ERROR
RMDcall = '''
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -J get_sra
#BSUB -q normal
#BSUB -outdir %s
#BSUB -o get_sra_%sJ.out
#BSUB -e get_sra_%sJ.screen
#BSUB -M 36000
#BSUB -n 1


echo "library(rmarkdown); rmarkdown::render('"%s"')" | R --no-save --no-restore

'''%(geo_dir,"%","%",pheno_file)


 
 if __name__ == '__main__':
 	#Make directory if does not exist
	make_dir(geo_dir)

	make_dir(data_dir)
    
	make_dir(res_dir)

	#Download matrix file in data directory if it does not exist
	#download_file(matrix_file_url, file_path)
	subprocess.call(mat_cmd, shell=True)
	print("Matrix file downloaded")

	#Write phenotype RMD file
	write_in_file(pheno_file,Rcall)
	print("Phenotype and analysis RMD file prepared")

	#Change demo to GEO ID in phenotype RMD file
	subprocess.call(cmd,shell=True)
	print("Demo changed to GEO ID entered in RMD file")

	#Change data and result dir paths in the phenotype RMD file
	subprocess.call(dat_cmd,shell=True)
	subprocess.call(res_cmd,shell=True)
	print("Data and Result directory path changes in RMD file")

	#Prepare lsf file for rmarkdown of phenotype file
	write_in_file("get_sra.lsf",RMDcall)
	print("get_sra file prepared")








