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
# source all your modules here
#source("modules/sidebar.R")


#######################################################################################################################
#######################################################################################################################
############################################       FUNCTIONS       ####################################################
#######################################################################################################################
#######################################################################################################################

# Sort file and tabix
sort_and_tabix <- function(tail_f) {
  bgzip <- file.path(HTSLIB_PATH,"bgzip")
  tabix <- file.path(HTSLIB_PATH,"tabix")
  
  tail_f_sorted <- gsub(tail_ext, ".sorted.csv", tail_f)
  tabix_l_file <- gsub(tail_ext, tabix_l_ext, tail_f)
  tail_f_gz <- paste0(tail_f_sorted, ".gz")
  sort_cmd=paste0("{ head -n 1 ",tail_f, " && tail -n +2 ", tail_f, " | sort -k",AGI_col," -k",start_col,"n,",end_col,"n; } > ", tail_f_sorted)
  tabix_cmd <- paste(bgzip, tail_f_sorted, "&&",  tabix, tail_f_gz, "-S 1 -s",AGI_col,"-b",start_col,"-e",end_col)
  
  tabix_l_cmd <- paste(tabix, "-l", tail_f_gz, ">", tabix_l_file)
  system(sort_cmd)
  system(tabix_cmd)
  
  system(tabix_l_cmd)
}

# Check if tabixed and tabix if not
check_if_tabixed <- function(tabix_file,tail_file, index=".tbi") {
  if (!file.exists(tail_file)) {
    shinyalert("OUCH", paste("THERE IS NO", tail_file), type = "error")
  }
  
  if (!file.exists(tabix_file) | !file.exists(paste0(tabix_file, index))){
    
    show_modal_spinner(spin = "fingerprint", text=paste("indexing", basename(file.path(tabix_file)),"\n this could take a while")) # show the modal window
    sort_and_tabix(tail_file)
    remove_modal_spinner() # remove it when done
    
  } 
}

# Read only part of files containing users' AGI 
read_tabixed_files <- function(file, AGI) {

  dt <- fread(cmd = paste(file.path(HTSLIB_PATH, "tabix"), file, AGI, "-h"), header = F)

  return(dt)
}

# Read GFF file for specified AGI
read_GFF_file <- function(gff, AGI) {
  dt <- fread(cmd = paste("grep", AGI, gff, "| grep -P 'mRNA|exon'"), col.names = gff_colnames) %>%
    mutate(orientation=case_when(strand=='-' ~ 0,
                                 strand=='+' ~ 1),
           feat_type =case_when(feature=="mRNA" ~ "gene",
                                feature=="exon" ~ "subgene"),
           ROI=paste0(seqnames, ":", start,"-", end))
  
  print(head(dt))
  return(dt)
}


getIntronsRetained <- function(intron_name, retained_introns){
  print(intron_name)
  print(retained_introns)
  intron_name %in% paste0("intron",unlist(strsplit(retained_introns, ":")))
  
}

build_coords_df <- function(AGI_df, GFF_DF, intron_cols) {
  
  if (!is.na(intron_cols)) {
    print(intron_cols)
    cols_to_keep <- c(intron_cols, "read_core_id", "mRNA", "retention_introns", "coords_in_read", "feat_id")
    
    AGI_df_coords <- AGI_df %>% dplyr::select(any_of(cols_to_keep)) %>%
      left_join(GFF_DF, by = c("mRNA"="AGI")) %>%
      separate(read_core_id, into=c("read_id", "chr", "read_start", "read_end"), sep = ',' , remove = F)
    
    # This data frame looks for the values of feat_id that are also column names (intron1 intron2 ...)
    # and checks whether it is retained or not.
    to_rem <- AGI_df_coords %>%
      filter(feat_id %in% colnames(.)) %>%
      rowwise() %>%
      mutate(retained = case_when(feature=="intron" ~ get(as.character(feat_id)),
                                  TRUE ~ FALSE)) %>%
      filter(retained==TRUE)
    
    AGI_df_coords<- AGI_df_coords %>%
      filter(feature!="intron") %>%
      mutate(retained=FALSE)
    
    transcripts_df <- AGI_df_coords %>%
      filter(feature=="mRNA") %>%
      mutate(feature="transcript",
             feat_type="transcript",
             feat_id="transcript1",
             start=as.numeric(read_start),
             end=as.numeric(read_end))
    
    
    AGI_df_coords <-  bind_rows(AGI_df_coords, to_rem, transcripts_df)%>%
      arrange(read_core_id, start, end) %>%
      dplyr::select(-c(intron_cols))
    
  } else {
    cols_to_keep <- c("read_core_id", "mRNA", "retention_introns", "coords_in_read", "feat_id")
    
    AGI_df_coords <- AGI_df %>% dplyr::select(any_of(cols_to_keep)) %>%
      left_join(GFF_DF, by = c("mRNA"="AGI")) %>%
      separate(read_core_id, into=c("read_id", "chr", "read_start", "read_end"), sep = ',' , remove = F)
    

    AGI_df_coords<- AGI_df_coords %>%
      filter(feature!="intron") %>%
      mutate(retained=FALSE)
    
    transcripts_df <- AGI_df_coords %>%
      filter(feature=="mRNA") %>%
      mutate(feature="transcript",
             feat_type="transcript",
             feat_id="transcript1",
             start=as.numeric(read_start),
             end=as.numeric(read_end))
    
    
    AGI_df_coords <-  bind_rows(AGI_df_coords, transcripts_df)%>%
      arrange(read_core_id, start, end) 
  }

  write_tsv(AGI_df_coords, "~/DATA/testAGI.tsv")

  return(AGI_df_coords)

}


#######################################################################################################################
#######################################################################################################################
############################################       GLOBALS       ######################################################
#######################################################################################################################
#######################################################################################################################

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