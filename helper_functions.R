#source("global.R")


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

#### SERVER ####


shinyInput <- function(FUN,id,num,...) {
  inputs <- character(num)
  for (i in seq_len(num)) {
    inputs[i] <- as.character(FUN(paste0(id,i),label=NULL,...))
  }
  inputs
}


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
  
  return(dt)
}


getIntronsRetained <- function(intron_name, retained_introns){
  intron_name %in% paste0("intron",unlist(strsplit(retained_introns, ":")))
  
}

build_coords_df <- function(AGI_df, GFF_DF, intron_cols) {

  AGI_df_coords <- AGI_df  %>%
    left_join(GFF_DF, by = c("mRNA"="AGI"), relationship = "many-to-many") %>%
    separate(read_core_id, into=c("read_id", "chr", "read_start", "read_end"), sep = ',' , remove = F, convert = T) %>%
    mutate(addtail_nchar =nchar(additional_tail))
  AGI_df_coords$orientation=unique(GFF_DF$orientation)
  

  transcripts_df <- AGI_df_coords %>%
    filter(feature=="mRNA") %>%
    mutate(feature="transcript",
           feat_type="transcript",
           feat_id="transcript1",
           start=as.numeric(read_start),
           end=as.numeric(read_end))
  
  if (!is.na(intron_cols)) {

    cols_to_keep <- c(intron_cols,"orientation", "read_core_id", "mRNA", "retention_introns", "coords_in_read",  "polya_length", "additional_tail", "origin", "feat_id", "feature", "feat_type")

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

    AGI_df_coords <- AGI_df_coords %>%dplyr::select(all_of(c(cols_to_keep)))%>%
      bind_rows(AGI_df_coords, to_rem, transcripts_df)%>%
      dplyr::select(-c(intron_cols))
    

    
  } else  if (is.na(intron_cols)) {
    cols_to_keep <- c("orientation", "read_core_id", "mRNA", "retention_introns", "coords_in_read", "polya_length", "additional_tail", "origin", "feat_id", "feature", "feat_type")

    AGI_df_coords <- AGI_df_coords 
    AGI_df_coords <-  bind_rows(AGI_df_coords, transcripts_df)


    AGI_df_coords<- AGI_df_coords %>%
      filter(feature!="intron") %>%
      mutate(retained=FALSE)

    AGI_df_coords <-  AGI_df_coords%>%
      dplyr::select(all_of(c(cols_to_keep)))%>%
      bind_rows(AGI_df_coords, transcripts_df)



      
  } else {
    stop("something went wrong with introns")
  }
  
  polya_df <- AGI_df_coords %>%
    filter(feature=="mRNA") %>%
    mutate(feature="polyA_tail",
           feat_type="tail",
           feat_id="polyA_tail1",
           start=case_when(orientation==1 ~ read_end+1,
                           orientation==0 ~ read_start-1),
           end= case_when(orientation==1 ~ read_end+ round(polya_length),
                          orientation==0 ~ read_start -round(polya_length)))
  
  addtail_df <- AGI_df_coords %>%
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

  
  AGI_df_coords <-  bind_rows(AGI_df_coords, polya_df, addtail_df)%>%
    arrange(read_core_id, start, end)
  
  write_tsv(AGI_df_coords, "~/DATA/testAGI.tsv")
  
  return(AGI_df_coords)
  
}

#### UI ####

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
