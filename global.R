# Global.R is loaded once at App start (not at browser refresh!)

# unload all libraries
#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

# load all libraries
library(shiny)
library(tidyverse)
library(shinydashboard)
library(data.table)
library(shinyFiles)
library(shinybusy)
library(shinyalert)
library(pbapply)
library(cowplot)
library(gggenes)
library(GenomicRanges)
library(shinyBS)
library(DT)

HTSLIB_PATH = "/shared/biotools/htslib/1.18/bin/"

datasets_path <- "/shared/home/jpeter/ShinyApps/RNAvis/data/FLEPseq_runs"
runs_infos <- read.table("/shared/home/jpeter/ShinyApps/RNAvis/data/Runs_infos.tsv", col.names =c("Run_name", "Run_infos"), sep = '\t', fill=T) 

flep_runs_df <- get_flepruns(datasets_path)
flep_runs <- setNames(flep_runs_df$flep_run, flep_runs_df$run_desc)

ref_gff <- "/shared/home/jpeter/Data/ReferenceGenomes/A_thaliana/Araport11/Araport11_GTF_genes_transposons.current.gtf"
gff_colnames <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

mapping_dir="2_Mapping"
tail_dir="4_Tail"
tail_ext=".read_info.result.merged.parts.csv"
mapping_ext=".bam.cov.txt"
mapping_cols=c("rname","startpos","endpos", "numreads", "covbases", "coverage", "meandepth","meanbaseq","meanmapq")
sample_table="barcode_correspondance.*$"

index_ext=".sorted.csv.gz"
tabix_l_ext=".list.tsv"

AGI_col=4
start_col=6
end_col=7
