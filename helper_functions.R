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
library(shinyscreenshot)
library(shinyhelper)
library(shinyWidgets)
library(pbapply)
library(cowplot)
library(GenomicRanges)
library(DT)
library(ggsci)
library(markdown)
library(cowplot)
library(rcartocolor)



get_flepruns <- function(datasets_path) {
  flep_runs <- list.dirs(datasets_path, recursive = F, full.names = T)
  flep_runs_df <- as.data.frame(flep_runs) %>%
    mutate(Run_bname = basename(flep_runs)) %>%
    left_join(runs_infos, by=c("Run_bname"= "Run_name")) %>%
    mutate(run_desc=
             case_when(
               !is.na(Run_infos) ~ paste(Run_bname, Run_infos, sep=":"),
               TRUE ~ Run_bname))
  return(flep_runs_df)
  
}


shinyInput <- function(FUN,id,num,...) {
  inputs <- character(num)
  for (i in seq_len(num)) {
    inputs[i] <- as.character(FUN(paste0(id,i),label=NULL,...))
  }
  inputs
}


# Read only part of files containing users' transcript 
read_tabixed_files_single_region <- function(file, transcript) {
  dt <- fread(cmd = paste(file.path(HTSLIB_PATH, "tabix"), file, transcript, "-h"), header = F)
  return(dt)
  }

# Read only part of files containing users' transcript 
read_tabixed_files_multiple_regions <- function(file, transcripts) {
  transcripts_string <- paste(unique(unlist(strsplit(transcripts, "\n"))), collapse = " ")

  dt <- fread(cmd = paste(file.path(HTSLIB_PATH, "tabix"), file, transcripts_string, "-h"), header = F)
  }

# Read GTF file for specified transcript
read_GTF_file <- function(gtf, transcript) {
  

  #dt <- fread(cmd = paste("grep", transcript, gtf, "| grep -P 'mRNA|exon|five_prime_UTR|three_prime_UTR'"), col.names = gtf_colnames) %>%
  dt <- fread(cmd = paste("grep", transcript, gtf, "| grep -P 'mRNA|exon'"), col.names = gtf_colnames)
  
  #print(transcript)
  if (nrow(dt)>0) {
    dt <- dt %>%
      mutate(transcript=transcript,
             orientation=case_when(strand=='-' ~ 0,
                                   strand=='+' ~ 1),
             feat_type =case_when(feature=="mRNA" ~ "gene",
                                  TRUE ~ "subgene"),
             ROI=paste0(seqnames, ":", start,"-", end)) %>%
      select(-c(attributes))
    
    return(dt)
    
  } else {
    print("not in GFF")
  }


  
}

getIntronsRetained <- function(intron_name, retained_introns){
  intron_name %in% paste0("intron",unlist(strsplit(retained_introns, ":")))
  
}

build_coords_df <- function(transcript_df, GTF_DF, intron_cols) {
  transcript_df_coords <- transcript_df  %>%
    left_join(GTF_DF, by = c("mRNA"="transcript"), relationship = "many-to-many") %>%
    mutate(addtail_nchar =nchar(additional_tail))
  transcript_df_coords$orientation=unique(GTF_DF$orientation)
  transcripts_df <- transcript_df_coords %>%
    filter(feature=="mRNA") %>%
    mutate(feature="transcript",
           feat_type="transcript",
           feat_id="transcript1",
           start=as.numeric(read_start),
           end=as.numeric(read_end))
  
  
  if (length(intron_cols)>0) {

    cols_to_keep <- c(intron_cols,"orientation", "read_core_id", "mRNA", "retention_introns", "coords_in_read",  "polya_length", "additional_tail", "origin", "feat_id", "feature", "feat_type", "sample_id")

    # This data frame looks for the values of feat_id that are also column names (intron1 intron2 ...)
    # and checks whether it is retained or not.
    to_rem <- transcript_df_coords %>%
      filter(feat_id %in% colnames(.)) %>%
      rowwise() %>%
      mutate(retained = case_when(feature=="intron" ~ get(as.character(feat_id)),
                                  TRUE ~ FALSE)) %>%
      filter(retained==TRUE)

    
    transcript_df_coords<- transcript_df_coords %>%
      filter(feature!="intron") %>%
      mutate(retained=FALSE)

    transcript_df_coords <- transcript_df_coords %>%dplyr::select(all_of(c(cols_to_keep)))%>%
      bind_rows(transcript_df_coords, to_rem, transcripts_df)%>%
      dplyr::select(-c(intron_cols))
    

    
  } else  if (length(intron_cols)==0) {
    cols_to_keep <- c("orientation", "read_core_id", "mRNA", "retention_introns", "coords_in_read", "polya_length", "additional_tail", "origin", "feat_id", "feature", "feat_type", "sample_id")

    transcript_df_coords <-  bind_rows(transcript_df_coords, transcripts_df)

    transcript_df_coords<- transcript_df_coords %>%
      filter(feature!="intron") %>%
      mutate(retained=FALSE)

    transcript_df_coords <-  transcript_df_coords%>%
      dplyr::select(all_of(c(cols_to_keep)))%>%
      bind_rows(transcript_df_coords, transcripts_df)

      
  } else {
    stop("something went wrong with introns")
  }
  
  polya_df <- transcript_df_coords %>%
    filter(feature=="mRNA") %>%
    mutate(feature="polyA_tail",
           feat_type="tail",
           feat_id="polyA_tail1",
           start=case_when(orientation==1 ~ read_end+1,
                           orientation==0 ~ read_start-1),
           end= case_when(orientation==1 ~ read_end+ round(polya_length),
                          orientation==0 ~ read_start -round(polya_length)))
  
  addtail_df <- transcript_df_coords %>%
    filter(feature=="mRNA") %>%
    mutate(feature="add_tail",
           feat_type="tail",
           feat_id="add_tail1",
           start=case_when(orientation==1 ~ read_end +round(polya_length),
                           orientation==0 ~ read_start -round(polya_length)),
           end=case_when(orientation==1 & !is.na(additional_tail) ~ read_end +round(polya_length)+addtail_nchar,
                         orientation==1 & is.na(additional_tail) ~ read_end + round(polya_length),
                         orientation==0 & !is.na(additional_tail) ~ read_start - round(polya_length) - addtail_nchar,
                         orientation==0 & is.na(additional_tail) ~ read_start - round(polya_length)))

  
  transcript_df_coords <-  bind_rows(transcript_df_coords, polya_df, addtail_df) 
  

  return(transcript_df_coords)
  
}


build_intronic_profile_for_plot <- function(filteredintron) {
  coords_df <- filteredintron 
  
  df_coords_gene <- coords_df %>% 
    filter(feat_type=="gene")
  genotypes <- unique(coords_df$origin)
  gene_name <- unique(df_coords_gene$mRNA)
  
  # # gene: this is to have the 'model' of the intron profile
  df_coords_gene <- df_coords_gene[!duplicated(df_coords_gene$mRNA, by=c("retention_introns", "origin")),]
  df_coords_gene <- df_coords_gene[rep(seq_len(nrow(df_coords_gene)), each = length(genotypes)), ]
  df_coords_gene$origin=genotypes
  # subgene (intron exon)
  df_coords_subgene <- coords_df %>% filter(feat_type=="subgene") %>%
    mutate(parent_start=unique(df_coords_gene$start),
           parent_stop=unique(df_coords_gene$end))
  
  # tail info
  df_coords_tails <- coords_df %>% filter(feat_type=="tail") %>%
    mutate(parent_start=unique(df_coords_gene$start),
           parent_stop=unique(df_coords_gene$end))
  
  
  # transcript
  df_coords_transcript <- coords_df %>% 
    filter(feat_type=="transcript") 
  
  
  # filter subgene to remove introns/exons that exist in annotation but not present in read (ie not sequenced)
  #df_coords_transcript_subgene <- df_coords_subgene
  df_coords_transcript_subgene <- df_coords_subgene%>%
    filter(end>=read_start,
           start<=read_end) %>%
    mutate(start=case_when((strand == '+' & start<read_start) ~ read_start,  # same instructions for each strand as when strand=='-', start is the smallest coord value but is biologically the end of the transcript
                           (strand == '-' & start<read_start) ~ read_start, 
                           TRUE ~start),
           end = case_when((strand == '+' & end> read_end) ~read_end,
                           (strand == '-' & end> read_end) ~read_end,
                           TRUE ~end))
  
  # combining dataframes
  df_coords_transcript_subgene_tail <- rbind(df_coords_transcript_subgene, df_coords_tails)
  df_coords_transcript_subgene_tail$parent_start <- NA
  df_coords_transcript_subgene_tail$parent_stop <- NA
  
  # adjusting gene size to allow subgene tail which is by def. not within gene boundaries
  df_coords_transcript_subgene_tail$parent_start <- min(df_coords_transcript_subgene_tail$start, na.rm = T) -1
  df_coords_transcript_subgene_tail$parent_stop <- max(df_coords_transcript_subgene_tail$end, na.rm=T) +1
  
  intronic_profile_for_plot <- list(df_coords_transcript = df_coords_transcript, 
                                    df_coords_transcript_subgene_tail = df_coords_transcript_subgene_tail, 
                                    df_coords_gene = df_coords_gene,
                                    df_coords_subgene=df_coords_subgene,
                                    gene_name=gene_name)
  
  #print(df_coords_gene)
  
  return(intronic_profile_for_plot)
  
} 

plot_polya_bulk <- function(dataset, plot_type) {
  dataset <- dataset%>%
    group_by(origin) %>%
    mutate(label=paste0("n = ", n()),
           origin=paste0(origin, " (", label, ")"))
  
  if (plot_type =="density") {
    polya_bulk <- ggplot(dataset, aes(x=polya_length, fill=origin)) +
      geom_density(alpha=0.5, color="black")+
      ggcustom_theme +
      theme(legend.position="bottom") +
      theme(legend.background = element_rect(linewidth=0.5, linetype="solid"))
    
  } else {
    polya_bulk <- ggplot(dataset, aes(x=polya_length, fill=origin, color=origin)) +
      geom_bar()+
      facet_wrap(~origin) + 
      ggcustom_theme +
      theme(legend.position="bottom") +
      theme(legend.background = element_rect(linewidth=0.5, linetype="solid")) 
  }
  
  polya_bulk <- polya_bulk + ggtitle("Poly(A) length distribution")

  return(polya_bulk)
}

plot_urid <- function(dataset, threshold) {
  urid_table <-dataset %>%
    mutate(uridylated=case_when(add_tail_pct_T>=threshold ~ "uridylated",
                                add_tail_pct_T<threshold ~ "not_uridylated",
                                is.na(add_tail_pct_T) ~ "no_additional_tail")) %>%
    group_by(origin, uridylated) %>%
    summarise(n=n())
  
  print(dataset %>% select(add_tail_pct_T))
  
  ggplot(urid_table, aes(x=origin, y=n, fill=origin))+
    geom_col(position='stack')+
    theme(legend.position="bottom")+
    facet_wrap(~uridylated) +
    ggcustom_theme +
    ggtitle("Number of reads by uridylation status") +
    theme(legend.position = "bottom")
}


check_transcripts_list <- function(user_gene_list, total_gene_list) {

  user_list <- fread(user_gene_list$datapath, header=F, col.names = "AGI") %>%
    filter(AGI %in% total_gene_list$AGI) %>%
    distinct(AGI) %>%
    pull(AGI)
  
  return(user_list)
  
  #correct_list <- 
  
}

#### UI ####

ggcustom_theme <- list(
  theme_bw(),
  scale_color_carto_d(),
  scale_fill_carto_d()
)

dropdownMenuCustom <-     function (..., type = c("messages", "notifications", "tasks"), 
                                    badgeStatus = "primary", icon = NULL, .list = NULL, customSentence = customSentence) 
{
  type <- match.arg(type)
  if (!is.null(badgeStatus)) shinydashboard:::validateStatus(badgeStatus)
  items <- c(list(...), .list)
  lapply(items, shinydashboard:::tagAssert, type = "li")
  dropdownClass <- paste0("dropdown ", type, "-menu")
  if (is.null(icon)) {
    icon <- switch(type, messages = shiny::icon("envelope"), 
                   notifications = shiny::icon("warning"), tasks = shiny::icon("tasks"))
  }
  numItems <- length(items)
  if (is.null(badgeStatus)) {
    badge <- NULL
  }
  else {
    badge <- tags$span(class = paste0("label label-", badgeStatus), 
                       numItems)
  }
  tags$li(
    class = dropdownClass, 
    a(
      href = "#", 
      class = "dropdown-toggle", 
      `data-toggle` = "dropdown", 
      icon, 
      badge
    ), 
    tags$ul(
      class = "dropdown-menu", 
      tags$li(
        class = "header", 
        customSentence(numItems, type)
      ), 
      tags$li(
        tags$ul(class = "menu", items)
      )
    )
  )
}

customSentence <- function(numItems, type) {
  paste("Feedback & suggestions")
}




