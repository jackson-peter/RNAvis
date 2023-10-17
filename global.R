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
runs_infos <- fread("/home/jpeter/ShinyApps/RNAvis/data/Runs_infos.tsv", header=F, col.names=c("Run_name", "Run_infos")) 

flep_runs <- list.dirs(datasets_path, recursive = F, full.names = T)
#names(flep_runs) <- lapply(flep_runs, basename)
flep_runs_df <- as.data.frame(flep_runs) %>%
  mutate(Run_bname = basename(flep_runs)) %>%
  left_join(runs_infos, by=c("Run_bname"= "Run_name")) %>%
  mutate(run_desc=
           case_when(
             !is.na(Run_infos) ~ paste(Run_bname, Run_infos, sep=":"),
             TRUE ~ Run_bname))
  

flep_runs <- setNames(flep_runs_df$flep_run, flep_runs_df$run_desc)



ref_gff <- "/home/jpeter/DATA/ReferenceGenomes/Athaliana/Araport11/Araport11_GTF_genes_transposons.Jul2023.gtf"
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


