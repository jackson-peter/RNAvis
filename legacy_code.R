# LEGACY CODE


# Sort file and tabix
sort_and_tabix <- function(tail_f) {
  # path to softs
  bgzip <- file.path(HTSLIB_PATH,"bgzip")
  tabix <- file.path(HTSLIB_PATH,"tabix")
  tail_f_sorted <- gsub(tail_ext, ".sorted.csv", tail_f)
  tabix_l_file <- gsub(tail_ext, tabix_l_ext, tail_f)
  tail_f_gz <- paste0(tail_f_sorted, ".gz")
  sort_cmd=paste0("{ head -n 1 ",tail_f, " && tail -n +2 ", tail_f, " | sort -k",transcript_col," -k",start_col,"n,",end_col,"n; } > ", tail_f_sorted)
  tabix_cmd <- paste(bgzip, "-I 9", tail_f_sorted, "&&",  tabix, tail_f_gz, "-S 1 -s",transcript_col,"-b",start_col,"-e",end_col)
  
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




output$sample_table2 = DT::renderDataTable({
  req(global$sample_corr$genotype)
  
  DT::datatable(cbind(Pick=shinyInput(checkboxInput,"srows_",length(global$sample_corr$genotype),value=TRUE,width=1), global$sample_corr),
                options = list(orderClasses = TRUE,
                               lengthMenu = c(5, 25, 50),
                               pageLength = 25 ,
                               drawCallback= JS(
                                 'function(settings) {
                                     Shiny.bindAll(this.api().table().node());}')
                ),selection='none',escape=F)
  
})


# ## Selected Genotypes
# genoSelect <- reactive({
#   rows=names(input)[grepl(pattern = "srows_",names(input))]
#   paste(unlist(lapply(rows,function(i){
#     if(input[[i]]==T){
#       index=as.numeric(substr(i,gregexpr(pattern = "_",i)[[1]]+1,nchar(i)))
#       
#       return(global$sample_corr$genotype[index])
#     }
#   })))
#   
# })
# 
