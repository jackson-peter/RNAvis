
message("beginning app")
source("helper_functions.R")
source("config.R")
source("global.R")



theme_set(theme_bw())

# UI ###########################################################
# header ----
header <- 
  dashboardHeader(title = HTML("RNAvis"),
                   disable = F,
                   dropdownMenuCustom( type = 'message',
                                       customSentence = "feedback & suggestions",
                                       messageItem(
                                         from = "jackson.peter",#'Feedback and suggestions',
                                         message =  "",#paste0("jackson.peter@ibmp-cnrs.unistra.fr" ),
                                         icon = icon("envelope"),
                                         href = "mailto:jackson.peter@ibmp-cnrs.unistra.fr"),
                                       icon = icon('comment'))
                   ) #/ dashboardheader

# sidebar ----
sidebar <-
  dashboardSidebar(
    hr(),
    sidebarMenu(id = "tabs",
                menuItem("Run Selection", tabName = 'runTab', icon = icon('gears'), selected=TRUE),
                menuItem("Transcript-Specific", tabName = "transcriptSpecificTab", icon=icon("minus")),
                menuItem("Transcripts List", tabName = "transcriptsListTab", icon=icon("list-ul")),
                menuItem("README", tabName = "readme", icon=icon("mortar-board"))

    ), #/ sidebarmenu
    hr(),
    screenshotButton()
  ) #/ dashboardsidebar


# body ----
body <-
  dashboardBody(
  tags$head(includeCSS("www/custom.css")),
  
  tabItems(
    # sidebar Tab: run selector ----
    tabItem(tabName = "runTab",
            # Inputs box ---
            box(width = 12, status = "primary", solidHeader = T, title="chose a FLEPseq run",
                fluidRow(column(selectizeInput("runSelection", inputId = 'runSelection', label=NULL, choices = c("Choose a run" = "", flep_runs), multiple = FALSE), width=9 ),
                         actionButton("SubmitRunSel", "Get the data!"), width=3)
                )%>% 
              helper(icon = "question",
                     colour = "grey",
                     type = "markdown",
                     content = "runSelection"),
            # /box
            # plots about run
            tabsetPanel(id = "run_tabsetPanel",
              tabPanel(value = "dl_tabPanel", h5("Download Sample Data"),
                       fluidRow(column(dataTableOutput("sample_table"), width=12)),
                       hr(),
                       em("Full dataset per sample is available for download."),
                       selectizeInput("download_sample_sel", inputId = 'download_sample_sel', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE)),
                       downloadButton("download_sample_data", "Download")
                       ),
              tabPanel(h5("Run mapping stats"),
                       h4("Number of reads by sample"),
                       dataTableOutput("smpl_reads"),
                       hr(style = "border-top: 1px solid #000000;"),
                       h4("Percentage of reads by sample"),
                       plotOutput("pctreads"),
                       hr(style = "border-top: 1px solid #000000;"),
                       h4("Number of detected genes by sample"),
                       plotOutput("numgenes"),
                       hr(style = "border-top: 1px solid #000000;")
                       ),
              tabPanel(h5("Poly(A) Distribution"),
                       #fluidRow(column(imageOutput("bulk_ig_global"), width=12))
                       
                       fluidRow(column(h4("Bulk distribution of Poly(A) tail lengths"), width=12),
                                column(imageOutput("bulk_polyA_global"), width=12, align="center"),
                                column(hr(style = "border-top: 1px solid #000000;"), width=12)
                                ),
                       
                       fluidRow(column(h4("Intergenic distribution of Poly(A) tail lengths"), width=12),
                                column(imageOutput("ig_polyA_global"), width=12, align="center"),
                                column(hr(style = "border-top: 1px solid #000000;"), width=12)
                       ),
                       
                       # fluidRow(column(h4("Cumulative distribution of Poly(A) tail lengths"), width=12),
                       #          column(imageOutput("cumul_polyA_global"), width=12, align="center"),
                       #          column(hr(style = "border-top: 1px solid #000000;"), width=12)
                       # ),
                       

                       ),
              # tabPanel(h5("cumulative Poly(A) Distribution"),
              #          imageOutput("cumul_polyA_global")
              #          )
              ) #/ tabsetpanel

    ), # /tabitem run selector
    
    
    # tabitem transcriptSpecificTab -----
    tabItem(tabName = "transcriptSpecificTab",
            box(width = 12, 
                status = "primary", 
                solidHeader = TRUE, 
                title="Select your Transcript of interest",
                column(width=9,
                       selectizeInput("transcript_sel", inputId = 'transcript_sel', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE))),
                column(width=3,
                       actionButton(inputId = "Submittranscript", label = "Get/Update Data")),
            )%>% 
              helper(icon = "question",
                     colour = "grey",
                     type = "markdown",
                     content = "transcript_sel"), 
            tabsetPanel(id="transctipt_tabsetpanel",
              tabPanel(h5("Transcript Overview"),
                       plotOutput("plot_gene"),
                       dataTableOutput("GTFtable_single")
                       ), # /tabpanel
              tabPanel(h5("Intronic Profiles"),
                       fluidRow(box(width=12, status="primary", solidHeader = T, title="Intronic profile selection", collapsible=T, collapsed=F,
                                    uiOutput("intron_profile_sel")
                                    ) # /box
                       ),
                       
                       plotOutput("plot_reads", width = "100%",
                                  dblclick = "plot_reads_dblclick",
                                  brush = brushOpts(
                                    id = "plot_reads_brush",
                                    resetOnNew = TRUE))%>% 
                         helper(icon = "question",
                                colour = "grey",
                                type = "markdown",
                                content = "plot_reads"),
                       em("You can zoom in on the plot by first selecting an area of the plot and then by double-clicking on it."),
                       em("Double-clicking on the plot with no selection resets the plot to original scale.")
                       
                       
                       ), # /tabpanel
              tabPanel(h5("FLEPseq Results"),
                       box(width = 12, status = "primary", solidHeader = TRUE, title="Select a subset of columns",
                           selectizeInput('singleTranscript_columnSel', "by default: all columns",choices=NULL, selected=NULL, multiple=T), # columns selector
                           actionButton(inputId = "getFlepTable_single",label = "Get FLEPseq2 table"),
                           downloadButton("download_FlepTable_single", "Download")
                           ),
                       dataTableOutput("FlepTable_single")

                       ), # /tabPanel
              tabPanel(h5("Poly(A) Tail"),
                       box(width=12, 
                           
                           plotOutput("polyaDistr_single",
                                  dblclick = "polyaDistr_single_dblclick",
                                  brush = brushOpts(
                                    id = "polyaDistr_single_brush",
                                    resetOnNew = TRUE)),
                           radioGroupButtons(
                             inputId = "polyaHist_or_polyaDistr_single",
                             label = "plot Type:",
                             choices = c(
                               `<i class='fa fa-area-chart'></i>` = "density",
                               `<i class='fa fa-bar-chart'></i>` = "bar"),
                               justified = TRUE,
                             selected = "density"),
                           em("You can zoom in on the plot by first selecting an area of the plot and then by double-clicking on it."),
                           em("Double-clicking on the plot with no selection resets the plot to original scale."),
                           ), #/ box
                       
                       ), # /tabpanel
              tabPanel(h5("Additional Tail"),
                       box(width=12,
                           sliderInput("urid_thresh_single", "Uridylation percentage threshold:", value = 60, min = 0, max = 100, width= "50%")%>%
                             helper(icon = "question",
                                    colour = "grey",
                                    type = "markdown",
                                    content = "urid_thresh"),
                           plotOutput("urid_single"))
                       )
              ) # / tabsetpanel

    ), #/ tabitem transcriptspecificTab
    
    # Tabbitem transcriptlistTab -----
    tabItem(tabName = "transcriptsListTab",    
            box(width = 12, 
                status = "primary", 
                solidHeader = TRUE, 
                title="Select your transcripts of interest",
                column(width=12,
                       selectizeInput("transcripts_sel", inputId = 'transcripts_sel', label = "select genes from this list", choices = NULL, selected = NULL, multiple = TRUE, options = list(create = FALSE))),
                column(width=12,
                       fileInput("transcripts_upload", label="or upload a list", multiple = F, buttonLabel = "Browse...", placeholder = "No file selected")),
                column(width=3,
                       actionButton(inputId = "Submittranscripts", label = "Get/Update Data")),
            )%>%
              helper(icon = "question",
                     colour = "grey",
                     type = "markdown",
                     content = "transcripts_sel"),
            tabsetPanel(id="transcripts_tabsetpanel",
              tabPanel(h5("Transcripts Overview"),
                       dataTableOutput("GTFtable_list")
              ), # /tabpanel
              tabPanel(h5("FLEPseq results"),
                       box(width = 12, status = "primary", solidHeader = TRUE, title="Select a subset of columns",
                           selectizeInput('transcriptList_columnSel', "by default: all columns",choices=NULL, selected=NULL, multiple=T), # columns selector
                           actionButton(inputId = "getFlepTable_list",label = "Get FLEPseq2 table"),
                           downloadButton("download_FlepTable_list", "Download")
                           ),
                       dataTableOutput("FlepTable_list")
              ), # / tabpanel
              tabPanel(h5("Poly(A) tail"),
                       box(width=12,
                           plotOutput(
                             "polyaDistr_list",
                             dblclick = "polyaDistr_list_dblclick",
                             brush = brushOpts(id = "polyaDistr_list_brush",
                                               resetOnNew = TRUE)
                         ),
                         radioGroupButtons(
                           inputId = "polyaHist_or_polyaDistr_list",
                           label = "Visualize:",
                           choices = c(
                             `<i class='fa fa-area-chart'></i>` = "density",
                             `<i class='fa fa-bar-chart'></i>` = "bar"
                           ),
                           justified = TRUE,
                           selected = "density"
                         ),
                         em("You can zoom in on the plot by first selecting an area of the plot and then by double-clicking on it."),
                         em("Double-clicking on the plot with no selection resets the plot to original scale."),
                         
                       ), #/ box
                       
              ), # /tabpanel
              tabPanel(h5("Additional Tail"),
                       box(width=12,
                           sliderInput("urid_thresh_list", "Uridylation percentage threshold:", value = 60, min = 0, max = 100, width= "50%")%>%
                             helper(icon = "question",
                                    colour = "grey",
                                    type = "markdown",
                                    content = "urid_thresh"),
                           plotOutput("urid_list")
                       )
              )
                       
              
            )
    
    ), #/ tabitem transcriptlistTab
    tabItem(tabName="readme",
            tabsetPanel(id="readme_tabsetpanel",
              fluidRow(includeMarkdown(README))
              )
            )
    

  ) #/ tabitems
) #/ dashboardbody 


# UI construction --------------------------------------------------------------
ui <- dashboardPage(skin="green",
  header, 
  sidebar,
  body
)

# SERVER ###########################################################
server <- function(input, output, session) {
  # This function stops the execution and fixes the bug requiring R termination when closing app
  session$onSessionEnded(function() {
    stopApp()
  })
  
  global <- reactiveValues(datapath = getwd())
  observe_helpers(withMathJax = TRUE, help_dir = DOC_DIR)
  ### REACTIVE DATA ------------------------------------------------------------
  
  # Single zoomable plot
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  ## Mapping data (coverage etc...)
  MAP_data <- eventReactive(input$SubmitRunSel,{
    map_files <- global$sample_corr$map_file
    message(map_files)
    names(map_files) <- global$sample_corr$genotype
    mapping_q <- rbindlist(lapply(map_files, fread, col.names=mapping_cols), idcol = "sample") 
    message(mapping_q)
    return(mapping_q)
  })
  
  ## gene list data 
  gene_list_data <- eventReactive(input$SubmitRunSel,{
    gene_list_f <- global$sample_corr$gene_list
    names(gene_list_f) <- global$sample_corr$genotype
    gene_list <- rbindlist(lapply(gene_list_f, fread, col.names="AGI"), idcol = "sample") 
    
    return(gene_list)
  })
  

  ## Dataset for list of transcripts
  transcripts_data <- eventReactive(input$Submittranscripts, {
    req(input$runSelection,
    isTruthy(input$transcripts_sel)|isTruthy(input$transcripts_upload))
    if (isTruthy(input$transcripts_sel)) {
      list_transcripts <- input$transcripts_sel
      
    } else if (isTruthy(input$transcripts_upload)) {
      print("list")
      list_transcripts <- check_transcripts_list(input$transcripts_upload, gene_list_data())
    }
    
    show_modal_spinner(spin = "fingerprint", text="fetching transcripts from list. This could take a while") # show the modal window
    GTF_DF <- rbindlist(lapply(list_transcripts, read_GTF_file, gtf=ref_gtf))
    
    tabixed_list <- global$sample_corr$tabix_file
    names(tabixed_list) <- global$sample_corr$genotype
    tabixed_df = setDT(as.list(tabixed_list))
    column_names <- names(fread(cmd = paste('zcat ',global$sample_corr$tabix_file[1],' | head -n 1')))
    column_names <- c('sample', column_names)
    transcripts_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files_multiple_regions, transcripts=list_transcripts), idcol = "sample")
    colnames(transcripts_DF) <- column_names
    #transcripts_DF <- transcripts_DF %>%
    #  filter(sample %in% genoSelect())
    updateSelectizeInput(session, "transcriptList_columnSel", choices=colnames(transcripts_DF))
    shinyalert("Nice!", paste(nrow(transcripts_DF), "transcripts in table"), type = "success")
    transcripts_DATA <- list(transcripts_DF = transcripts_DF, GTF_DF = GTF_DF)
    
    remove_modal_spinner() # remove modal spinner when done
    return(transcripts_DATA)

  })
  
  ## transcript specific data (FLEPseq results, gtf and coordinates)
  transcript_data <- eventReactive(input$Submittranscript,{
    req(input$runSelection, input$transcript_sel)
    # GTF gymnastics to handle introns
    show_modal_spinner(spin = "fingerprint", text="fetching transcripts from list\n this could take a while") # show the modal window
    gtf_infos <- read_GTF_file(ref_gtf, input$transcript_sel)
    mRNA_df <- gtf_infos %>% filter(feature=="mRNA")
    exon_df <- gtf_infos %>% filter(feature=="exon")
    mRNA_GR <- GRanges(mRNA_df$ROI)
    exon_GR <- GRanges(exon_df$ROI)
    # if several exons(and therefore introns)
    if (nrow(exon_df)>1) {
      introns_regions <- as.data.frame(GenomicRanges::setdiff(mRNA_GR, exon_GR, ignore.strand=TRUE)) %>%
        mutate(feature="intron",
               feat_type="subgene",
               ROI=paste0(seqnames, ":", start,"-", end),
               orientation=unique(gtf_infos$orientation),
               strand=unique(gtf_infos$strand),
               source=unique(gtf_infos$source))
      
      GTF_DF <- rbind(gtf_infos, introns_regions, fill=T) %>%
        group_by(feature) %>%
        mutate(transcript=input$transcript_sel,
               feat_id=case_when(orientation==0 ~paste0(feature, rev(seq_along(feature))),
                                 orientation==1 ~paste0(feature, seq_along(feature))))%>%
        arrange(start)
    # if only one exon
    } else {
      GTF_DF <- gtf_infos %>%
        mutate(transcript=input$transcript_sel,
               feat_id=paste0(feature,"1")) %>%
        arrange(start)
    }
    
    # Building transcript dataset
    column_names <- names(fread(cmd = paste('zcat ',global$sample_corr$tabix_file[1],' | head -n 1')))
    column_names <- c('sample', column_names)
    tabixed_list <- global$sample_corr$tabix_file
    names(tabixed_list) <- global$sample_corr$genotype
    tabixed_df = setDT(as.list(tabixed_list))
    tryCatch(
      {
        transcript_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files_single_region, transcript=input$transcript_sel), idcol = "sample")
        shinyalert("Nice!", paste(nrow(transcript_DF), "transcripts in table"), type = "success")
      },
      error = function(cond) {shinyalert("OUCH", paste("Error:", cond), type = "error")}
    )
    
    tryCatch(
      expr = {
        colnames(transcript_DF) <- column_names
        print(column_names)
        gff_strand=unique(GTF_DF$strand)
        transcript_DF <- transcript_DF %>%
          separate(read_core_id, into=c("read_id", "chr", "read_start", "read_end"), sep = ',' , remove = F, convert = T)%>%
          mutate(gff_strand=gff_strand) %>%
          mutate(
                 sort_val1 =case_when(
                   gff_strand=='+'~ read_start,
                   gff_strand=='-' ~ read_end),
                 sort_val2 =case_when(
                   gff_strand=='+' ~ read_end+(round(polya_length) + nchar(additional_tail)),
                   gff_strand=='-' ~ read_start - (round(polya_length) + nchar(additional_tail)),
                 ))%>%
          arrange(sample, sort_val1, sort_val2) %>%
          select(-c(sort_val1, sort_val2, gff_strand)) %>%
          group_by(sample, retention_introns) %>% 
          mutate(sample_id = row_number()) %>% 
          ungroup() %>%
          mutate(Run=basename(global$datapath))%>%
          mutate(across(retention_introns, as.character))%>%
          mutate(retention_introns = replace_na(retention_introns, "none")) 
        
        setDT(transcript_DF)
        

        #%>% filter(sample %in% genoSelect())
        
        n_introns <- unique(transcript_DF$mRNA_intron_num)
        if (n_introns>0) {
          intron_cols <- paste0("intron",1:n_introns)
          
          transcript_DF <- transcript_DF[, paste(intron_cols) := lapply(paste(intron_cols), getIntronsRetained, retained_introns = as.character(retention_introns)), by = retention_introns]
          coords_df <- build_coords_df(transcript_DF, GTF_DF, intron_cols) 
        } else {
          coords_df <- build_coords_df(transcript_DF, GTF_DF, vector()) 
        }

        updateSelectizeInput(session, "singleTranscript_columnSel", choices=colnames(transcript_DF))
      },
      error = function(e){
        shinyalert("Erf!", paste("Something went wrong", e), type = "error")
        print(e)}
     
    )


    transcript_DATA <- list(transcript_DF = transcript_DF, GTF_DF = GTF_DF, COORDS_DF = coords_df)
    remove_modal_spinner() # remove modal spinner when done

    return(transcript_DATA)
  })

  ## Dataset filtered by intronic profile
  filteredintron <- reactive({ 
    
    filtered_df<- transcript_data()$COORDS_DF %>% 
        dplyr::filter(retention_introns %in% as.vector(input$retention_introns)) 
    
    return(filtered_df)


  })
  
  ## Dataset column subsets (1 transcript) ----
  filtereddata_single <- eventReactive(input$getFlepTable_single,{
    
    if (is.null(input$singleTranscript_columnSel)) {
      filtereddata <- transcript_data()$transcript_DF
    } else {
      filtereddata <- transcript_data()$transcript_DF%>% 
        ungroup()%>%
        select(input$singleTranscript_columnSel)
    }
    return(filtereddata)
    
  })
  
  ## Dataset column subsets (multiple transcripts) ----
  filtereddata_list <- eventReactive(input$getFlepTable_list,{
    
    if (is.null(input$transcriptList_columnSel)) {
      filtereddata <- transcripts_data()$transcripts_DF
    } else {
      filtereddata <- transcripts_data()$transcripts_DF%>% 
        ungroup()%>%
        select(input$transcriptList_columnSel)
    }
    return(filtereddata)
    
  })
  
  ### OBSERVERS ---------------------------------------------------------------

  ## Observer for runSelection. 
  observeEvent(ignoreNULL = TRUE,
               #eventExpr = {input$runSelection},
               eventExpr = {input$SubmitRunSel},
               handlerExpr = {
                 
                 # Sets paths variables
                 
                 if (!input$runSelection == "") {

                   global$datapath <- input$runSelection
                   global$bulk_polyA_global_f <- file.path(global$datapath, paste0(basename(global$datapath), polya_bulk_plot_ext))
                   global$ig_polyA_global_f <- file.path(global$datapath, paste0(basename(global$datapath), polya_ig_plot_ext))
                   global$bulk_ig_global_f <- file.path(global$datapath, paste0(basename(global$datapath), bulk_ig_plot_ext))
                   global$cumul_polyA_global_f <- file.path(global$datapath, paste0(basename(global$datapath), cumul_polyA_plot_ext))
                   tabix_files <- list.files(path=file.path(global$datapath,tail_dir), pattern = paste0(index_ext, "$"),full.names = T)

                   sample_file <- list.files(path=global$datapath, pattern=sample_table, full.names = T)

                   mapping_files <- list.files(path=file.path(global$datapath,mapping_dir), pattern=paste0(mapping_ext, "$"), full.names=T)
                   print(sample_file)
                   global$sample_corr <- fread(sample_file, col.names = c("sample", "genotype"), header = F) %>%
                     mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(sample, index_ext)),
                            map_file=file.path(global$datapath, mapping_dir, paste0(sample, mapping_ext)),
                            gene_list=file.path(global$datapath, tail_dir, paste0(sample, tabix_l_ext))
                            )
                   print("ok")
                   
                     

                   tabix_files_txt <- paste(basename(global$sample_corr$tabix_file), collapse='\n')
                   shinyalert("Nice!", paste("Successfully added", tabix_files_txt, sep="\n"), type = "success")
                   genes_list=unique(rbindlist(lapply(global$sample_corr$gene_list, fread, header=F)))
                   updateSelectizeInput(session, 'transcript_sel', label=NULL, selected="", choices = genes_list$V1, options = list(create = FALSE), server = TRUE)
                   updateSelectizeInput(session, 'transcripts_sel', label=NULL, selected="", choices = genes_list$V1, options = list(create = FALSE), server = TRUE)
                   listf <- global$sample_corr$tabix_file
                   names(listf) <- basename(global$sample_corr$tabix_file)
                   updateSelectizeInput(session, 'download_sample_sel', label=NULL, selected="", choices = listf, options = list(create = FALSE), server = TRUE)
                   updateTabsetPanel(session, "run_tabsetPanel",
                                     selected = "dl_tabPanel")
                 }
               }) # end observer runSelection
  
  

  ### OUTPUTS ---------------------------------------------------------------

  ## intronic profile selection (radiobuttons) ----
  output$intron_profile_sel <- renderUI({
    req(transcript_data())
    tryCatch( { 
      selectizeInput('retention_introns', 'Select Retention intron', choices =unique(transcript_data()$transcript_DF$retention_introns), selected="none")
    
    }, error=function(e) {
      selectizeInput('retention_introns', 'Select Retention intron', choices =unique(transcript_data()$transcript_DF$retention_introns))
    })
  })
  
  output$sample_table2 = DT::renderDataTable({
    req(global$sample_corr$genotype)
    
    DT::datatable(global$sample_corr)
    
  })
  
  output$sample_table = DT::renderDataTable({
    req(global$sample_corr$genotype)
    out_table <- global$sample_corr %>%
      mutate(tabix_file= basename(tabix_file),
             map_file=basename(map_file),
             gene_list=basename(gene_list))
    
    DT::datatable(out_table)
    
  })
  
  output$download_sample_data <- downloadHandler(
    filename <- function() {
      basename(input$download_sample_sel)
    },
    
    content <- function(file) {
      file.copy(input$download_sample_sel, file)
    },
    contentType = NULL
  )
  
  
  
  ## table with user-selected column (1 transcript) ----
  output$FlepTable_single  <- renderDataTable(filtereddata_single(),options = list(scrollX = T))
  
  output$download_FlepTable_single <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0(basename(input$runSelection),"_", input$transcript_sel, ".tsv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write_tsv(filtereddata_single(), file)
    }
  )
  
  ## table with user-selected column (multiple transcripts) ----
  output$FlepTable_list  <- renderDataTable(filtereddata_list(),options = list(scrollX = T))
  
  output$download_FlepTable_list <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0(basename(input$runSelection),"_multiple_transcripts.tsv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write_tsv(filtereddata_list(), file)
    }
  )
  
  ## GTF table ----
  output$GTFtable_single <- renderDataTable(
    req(transcript_data()$GTF_DF %>%
          relocate(transcript)),
    options = list(pageLength = 50)
  )
  
  output$GTFtable_list <- renderDataTable(
    req(transcripts_data()$GTF_DF %>%
          relocate(transcript)),
    options = list(pageLength = 50)
  )
  
  ## coverage plot ----
  output$coverage <- renderPlot({
    req(input$SubmitRunSel)
    print(MAP_data())
    cov <- MAP_data() %>%
      group_by(sample, rname) %>%
      summarise(mean_cov=mean(coverage),
                mean_mapq =mean(meanmapq))
    ggplot(cov, aes(x=rname, y=mean_cov, fill=sample)) +
      geom_col(position = "dodge", alpha=0.5, color="black") +
      ggtitle("Mean coverage") +
      ggcustom_theme +
      theme(
        legend.position.inside = c(.05, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")
      ) +
      xlab("Chromosome") +
      ylab("Mean Coverage")
  })

  ## nreads plot ----
  output$numgenes <- renderPlot({
    req(input$SubmitRunSel)
    print(gene_list_data())

    genes_by_sample <- gene_list_data() %>%
      group_by(sample) %>%
      dplyr::summarise(n_genes=n())
    

    ggplot(genes_by_sample, aes(x=sample, y=n_genes, fill=sample)) +
      geom_col(position = "dodge", alpha=0.5, color="black") +
        #ggtitle("Number of genes by sample") +
        ggcustom_theme +
        theme(
          legend.position = "none",
          axis.text.x=element_text(angle=90, hjust=1)) +
        xlab("sample") +
        ylab("Number of genes")

  })
  
  
  output$smpl_reads <- DT::renderDataTable({
    req(input$SubmitRunSel)
    
    cov <- DT::datatable(MAP_data() %>%
      select(-c("startpos", "endpos", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")) %>%
      group_by(sample)%>%
      mutate(tot_read_smpl=sum(numreads)) %>%
      print() %>%
      pivot_wider(names_from = rname, values_from = numreads) %>%
      
      dplyr::rename(
        "Sample" = sample,
        "Number of reads"= tot_read_smpl)) 
  })
  
  
  ## pctreads plot ----
  output$pctreads <- renderPlot({
    req(input$SubmitRunSel)
    
    cov <- MAP_data() %>%
      group_by(sample)%>%
      mutate(tot_read_smpl=sum(numreads)) %>%
      group_by(sample, rname) %>%
      summarise(pct_chr=numreads/tot_read_smpl,
                mean_nreads=mean(numreads),
                mean_mapq =mean(meanmapq))
    ggplot(cov, aes(x=rname, y=pct_chr, fill=sample)) +
      geom_col(position = "dodge", alpha=0.5, color="black") +
      #ggtitle("Percentage of reads by sample") +
      ggcustom_theme +
      theme(
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")
      ) +
      xlab("Chromosome") +
      ylab("Read Percentage") +
      scale_y_continuous(labels = scales::percent)
  })
  
  
  
  ## mapping quality plot ----
  output$mapq <- renderPlot({
    req(input$SubmitRunSel)
    cov <- MAP_data() %>%
      group_by(sample, rname) %>%
      summarise(mean_cov=mean(coverage),
                mean_mapq =mean(meanmapq))
    ggplot(cov, aes(x=rname, y=mean_mapq, fill=sample)) +
      geom_col(position = "dodge", alpha=0.5, color="black") +
      ggtitle("Mean Mapping Quality")+
      ggcustom_theme +
      theme(legend.position = "none") +
      xlab("Chromosome") +
      ylab("Mean Mapping Quality")
  })
  
  ## polyA bulk distribution for single transcript ----
  output$polyaDistr_single <- renderPlot({
    req(transcript_data())
    polya_bulk <- plot_polya_bulk(transcript_data()$transcript_DF, input$polyaHist_or_polyaDistr_single)
    polya_bulk+
      coord_cartesian(xlim = ranges$x_polyaDistr_single, ylim = ranges$y_polyaDistr_single, expand = FALSE)
    

  })
  
  ## # uridylated reads for one transcript ----
  output$urid_single <- renderPlot({
    req(transcript_data())
    urid <- plot_urid(transcript_data()$transcript_DF, input$urid_thresh_single)
    urid
    
  })
  
  ## # uridylated reads for multiple transcript ----
  output$urid_list <- renderPlot({
    req(transcripts_data())
    urid <- plot_urid(transcripts_data()$transcripts_DF, input$urid_thresh_list)
    urid
    
  })
  
  ## polyA bulk distribution for multiple transcripts ----
  output$polyaDistr_list <- renderPlot({
    req(transcripts_data())
    polya_bulk <- plot_polya_bulk(transcripts_data()$transcripts_DF, input$polyaHist_or_polyaDistr_list)
    polya_bulk+
      coord_cartesian(xlim = ranges$x_polyaDistr_list, ylim = ranges$y_polyaDistr_list, expand = FALSE)
    
  })
  
  ## gene plot ----
  output$plot_gene <- renderPlot({
    req(transcript_data())
    GTF_DF <- transcript_data()$GTF_DF
    gene_size <- max(abs(GTF_DF$start - GTF_DF$end))
    chromo <- unique(GTF_DF$seqnames)
    gene <- unique(GTF_DF$transcript)
    strand <- unique(GTF_DF$strand)
    GTF_DF <- GTF_DF%>%
      filter(feature!="mRNA") %>%
      mutate(feat_id2= case_when(startsWith(feat_id, "exon") ~ "",
                                TRUE ~ gsub("\\D", "", feat_id)))
    
    y_max <- round(gene_size/50)
    
    ggplot(GTF_DF) +
      geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=y_max, fill=feature)) +
      geom_text(aes(x=start+(end-start)/2, y=0+y_max/2, label=feat_id2), size=3, color="white") +
      coord_fixed() +
      ggcustom_theme +
      theme(axis.text.y=element_blank(), 
            axis.ticks.y=element_blank(),
            legend.position="bottom",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab(chromo) + ylab("") +
      ggtitle(paste(gene, "(", strand, "strand)")) 
    
  })

  output$bulk_polyA_global <- renderImage({
    req(input$SubmitRunSel)
    # Return a list containing the filename
    list(src = global$bulk_polyA_global_f, style="height: 100%")
    }, deleteFile = FALSE)
  
  output$ig_polyA_global <- renderImage({
    req(input$SubmitRunSel)
    list(src=global$ig_polyA_global_f, style="height: 100%")
  }, deleteFile = FALSE)
  
  output$bulk_ig_global <- renderImage({
    req(input$SubmitRunSel)
    #list(src=global$bulk_ig_global_f, width = "100%", height = "700")
    list(src=global$bulk_ig_global_f, style="height: 100%")
  }, deleteFile = FALSE)
  
  output$cumul_polyA_global <- renderImage({
    req(input$SubmitRunSel)
    # Return a list containing the filename
    list(src = global$cumul_polyA_global_f, style="height: 100%")
  }, deleteFile = FALSE)

  
  
  ## plot reads (intronic profile) ----
  output$plot_reads <- renderPlot({
    req(transcript_data())
    
    validate(need(!is.null(input$retention_introns), "Please select an intron profile"))
    total_data <- build_intronic_profile_for_plot(filteredintron())

    # plot
    ggplot() +
      # transcripts.....
      #geom_rect(data=total_data$df_coords_transcript,
      #          aes(xmin = start, xmax = end, ymin = sample_id, ymax=sample_id + 0.5 ), fill="grey",alpha=.8)+
      geom_rect(data = total_data$df_coords_transcript_subgene_tail,
                aes(xmin = start, xmax = end, ymin = sample_id , ymax=sample_id +0.5,fill = feature), alpha=.8 ) +
      geom_segment(data=total_data$df_coords_transcript, aes(x = start, y = sample_id+0.25, xend = end, yend = sample_id+0.25),color = "black") +
      # .... annotated
      
      
      # model gene ...
      # geom_rect(data = total_data$df_coords_subgene,
      #           aes(xmin = start, xmax = end, ymin = 0, ymax= 0.5 , fill = feature)) +
      # geom_segment(data=total_data$df_coords_gene, 
      #                 aes(x = start, xend = end, y = 0.25, yend= 0.25), color = "black") +
      # # ... annotated
      
       
      facet_wrap(~sample, ncol = 1) +
      theme(strip.background =element_rect(fill="darkgrey"))+
      theme(strip.text = element_text(colour = 'white')) +
      theme(legend.position="bottom") +
      theme(legend.background = element_rect(linewidth=0.5, linetype="solid")) +
      ggtitle(total_data$gene_name) +
      coord_cartesian(xlim = ranges$x_plot_reads, ylim = ranges$y_plot_reads, expand = T) +
      theme(strip.text = element_text(
        size = 20, color = "black"))
    
  }, height = 1500)
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$polyaDistr_list_dblclick, {
    polyaDistr_list_brush <- input$polyaDistr_list_brush
    if (!is.null(polyaDistr_list_brush)) {
      ranges$x_polyaDistr_list <- c(polyaDistr_list_brush$xmin, polyaDistr_list_brush$xmax)
      ranges$y_polyaDistr_list <- c(polyaDistr_list_brush$ymin, polyaDistr_list_brush$ymax)
      
    } else {
      ranges$x_polyaDistr_list <- NULL
      ranges$y_polyaDistr_list <- NULL
    }
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot_reads_dblclick, {
    print("GOT IT")
    plot_reads_brush <- input$plot_reads_brush
    if (!is.null(plot_reads_brush)) {
      ranges$x_plot_reads <- c(plot_reads_brush$xmin, plot_reads_brush$xmax)
      ranges$y_plot_reads <- c(plot_reads_brush$ymin, plot_reads_brush$ymax)
      
    } else {
      ranges$x_plot_reads <- NULL
      ranges$y_plot_reads <- NULL
    }
    print(ranges$x_plot_reads)
    print(ranges$y_plot_reads)
    return(ranges)
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$polyaDistr_single_dblclick, {
    polyaDistr_single_brush <- input$polyaDistr_single_brush
    if (!is.null(polyaDistr_single_brush)) {
      ranges$x_polyaDistr_single <- c(polyaDistr_single_brush$xmin, polyaDistr_single_brush$xmax)
      ranges$y_polyaDistr_single <- c(polyaDistr_single_brush$ymin, polyaDistr_single_brush$ymax)
      
    } else {
      ranges$x_polyaDistr_single <- NULL
      ranges$y_polyaDistr_single <- NULL
    }
  })
 
} # /server

# App calling ------------
shinyApp(ui, server)

