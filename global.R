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

HTSLIB_PATH = "/biotools/htslib/1.9/bin/"

datasets_path <- "/home/jpeter/ShinyApps/RNAvis/data"
flep_runs <- list.dirs(datasets_path, recursive = F, full.names = T)
names(flep_runs) <- lapply(flep_runs, basename)
runs_infos <- fread("/home/jpeter/ShinyApps/RNAvis/data/Runs_infos.tsv", header=F, col.names=c("Run_name", "Run_infos")) 

ref_gff <- "/home/jpeter/DATA/ReferenceGenomes/Athaliana/Araport11/Araport11_GTF_genes_transposons.Jul2023.gtf"
gff_colnames <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

mapping_dir="2_Mapping"
tail_dir="4_Tail"
tail_ext=".read_info.result.merged.parts.csv"
mapping_ext=".bam.cov.txt"
mapping_cols=c("rname","startpos","endpos", "numreads", "covbases", "coverage", "meandepth","meanbaseq","meanmapq")
sample_table="barcode_correspondance.t.*"

index_ext=".sorted.csv.gz"
tabix_l_ext=".list.tsv"

AGI_col=4
start_col=6
end_col=7


