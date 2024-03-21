source("helper_functions.R")
source("global.R")



theme_set(theme_bw())
print(README)

# UI ###########################################################
# header =======================================================
header <- 
  dashboardHeader( title = HTML("RNAvis"),
                   disable = F,
                   dropdownMenuCustom( type = 'message',
                                       customSentence = "feedback & suggestions",
                                       # feedback notif
                                       messageItem(
                                         from = "jackson.peter",#'Feedback and suggestions',
                                         message =  "",#paste0("jackson.peter@ibmp-cnrs.unistra.fr" ),
                                         icon = icon("envelope"),
                                         href = "mailto:jackson.peter@ibmp-cnrs.unistra.fr"),
                                       icon = icon('comment'))
                   ) #/ dashboardheader

# sidebar =======================================================
sidebar <-
  dashboardSidebar(
    hr(),
    sidebarMenu(id = "tabs",
                menuItem("Run Selection", tabName = 'runTab', icon = icon('gears'), selected=TRUE),
                menuItem("Transcript-Specific", tabName = "transcriptSpecificTab", icon=icon("minus")),
                menuItem("Transcripts List", tabName = "transcriptsListTab", icon=icon("list-ul")),
                menuItem("README", tabName = "readme", icon=icon("mortar-board"))
                
                #menuItem("About", tabName = "about", icon = icon("question"))
                
    ), #/ sidebarmenu
    hr(),
    screenshotButton()
  ) #/ dashboardsidebar


# body =======================================================
body <-
  dashboardBody(
  tags$head(includeCSS("www/custom.css")),
  # Sidebar Tab Main dashboard -----
  tabItems(
    # sidebar Tab: run selector
    tabItem(tabName = "runTab",
            # Inputs box ---
            box(width = 12, status = "primary", solidHeader = TRUE, title="chose a FLEPseq run",
                selectizeInput("runSelection", inputId = 'runSelection', label=NULL, choices = c("Choose a run" = "", flep_runs), multiple = FALSE),
                fluidRow(column(dataTableOutput("sample_table"), width=11)),
                fluidRow(column(textOutput("dirname"),width=11)),
                fluidRow(column(textOutput("genotypes"),width=11))
                ),# /box
            # plots about run
            box(width = 12, status = "primary", solidHeader = TRUE, title="Run mapping & Coverage", collapsible = T, collapsed=T,
                fluidRow(column(splitLayout(plotOutput("coverage"),
                                            plotOutput("mapq"), cellWidths = c("50%", "50%"))
                                ,width=11))
                ) #/ box

    ), # /tabitem dashboard
    
    # Sidebar tab transcriptSpecificTab -----
    tabItem(tabName = "transcriptSpecificTab",
            box(width = 12, 
                status = "primary", 
                solidHeader = TRUE, 
                title="Select your Transcript of interest",
                column(width=12,
                       selectizeInput("transcript_sel", inputId = 'transcript_sel', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE))),
                column(width=3,
                       actionButton(inputId = "Submittranscript", label = "Get Data")),
            ), # /box
            tabsetPanel(
              tabPanel(h5("Transcript Overview"),
                       plotOutput("plot_gene"),
                       dataTableOutput("GTFtable_single")
                       ), # /tabpanel
              tabPanel(h5("Intronic Profiles"),
                       fluidRow(box(width=12, status="primary", solidHeader = T, title="Intronic profile selection", collapsible=T, collapsed=F,
                                    uiOutput("checkbox_introns")
                                    ) # /box
                       ),
                       em("You can zoom in on the plot by first selecting an area of the plot and then by double-clicking on it."),
                       plotOutput("plot_reads", width = "100%",
                                  dblclick = "plot_reads_dblclick",
                                  brush = brushOpts(
                                    id = "plot_reads_brush",
                                    resetOnNew = TRUE)
                                  )
                       
                       
                       ), # /tabpanel
              tabPanel(h5("FLEPseq Results"),
                       box(width = 12, status = "primary", solidHeader = TRUE, title="Select a subset of columns",
                           selectizeInput('singleTranscript_columnSel', "by default: all columns",choices=NULL, selected=NULL, multiple=T), # columns selector
                           actionButton(inputId = "getFlepTable_single",label = "Get FLEPseq2 table"),
                           downloadButton("download_FlepTable_single", "Download")
                           ),
                       dataTableOutput("FlepTable_single")

                       ), # /tabPanel
              tabPanel(h5("PolyA Distribution"),
                       plotOutput("polyA_distr_transcript")
                       ) # /tabpanel
              ) # / tabsetpanel

    ), #/ tabitem transcriptspecificTab
    
    # Sidebar tab transcriptlistTab -----
    tabItem(tabName = "transcriptsListTab",    
            box(width = 12, 
                status = "primary", 
                solidHeader = TRUE, 
                title="Select your transcripts of interest",
                column(width=12,
                       selectizeInput("transcripts_sel", inputId = 'transcripts_sel', label = NULL, choices = NULL, selected = NULL, multiple = TRUE, options = list(create = FALSE))),
                column(width=3,
                       actionButton(inputId = "Submittranscripts", label = "Get Data")),
            ), # /box
            tabsetPanel(
              tabPanel(h5("Transcripts Overview"),
                       #textOutput("test1"),
                       dataTableOutput("GTFtable_list")
              ), # /tabpanel
              tabPanel(h5("FLEPseq results"),
                       box(width = 12, status = "primary", solidHeader = TRUE, title="Select a subset of columns",
                           selectizeInput('transcriptList_columnSel', "by default: all columns",choices=NULL, selected=NULL, multiple=T), # columns selector
                           actionButton(inputId = "getFlepTable_list",label = "Get FLEPseq2 table"),
                           downloadButton("download_FlepTable_list", "Download")
                           ),
                       dataTableOutput("FlepTable_list"),
              ), # / tabpanel
              tabPanel(h5("PolyA Distribution"),
                       plotOutput("polyA_distr_transcripts")
              ) # /tabpanel
            )
    
    ), #/ tabitem transcriptlistTab
    tabItem(tabName="readme",
            tabsetPanel(
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
  
  ### REACTIVE DATA ------------------------------------------------------------
  
  # Single zoomable plot
  ranges <- reactiveValues(x = NULL, y = NULL)

  ## Mapping data (coverage etc...)
  MAP_data <- eventReactive(input$runSelection,{
    req(global$sample_corr)
    map_files <- global$sample_corr$map_file
    names(map_files) <- global$sample_corr$genotype
    mapping_q <- rbindlist(lapply(map_files, fread, col.names=mapping_cols), idcol = "origin") %>%
      filter(origin %in% genoSelect())
  })
  
  ## Selected Genotypes
  genoSelect <- reactive({
    rows=names(input)[grepl(pattern = "srows_",names(input))]
    paste(unlist(lapply(rows,function(i){
      if(input[[i]]==T){
        index=as.numeric(substr(i,gregexpr(pattern = "_",i)[[1]]+1,nchar(i)))
        
        return(global$sample_corr$genotype[index])
      }
    })))
    
  })
  
  ## Dataset for list of transcripts
  transcripts_data <- eventReactive(input$Submittranscripts, {
    req(input$runSelection, input$transcripts_sel)
    GTF_DF <- rbindlist(lapply(input$transcripts_sel, read_GTF_file, gtf=ref_gtf))
    tabixed_list <- global$sample_corr$tabix_file
    names(tabixed_list) <- global$sample_corr$genotype
    tabixed_df = setDT(as.list(tabixed_list))
    column_names <- names(fread(cmd = paste('zcat ',global$sample_corr$tabix_file[1],' | head -n 1')))
    column_names <- c('origin', column_names)
    transcripts_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files_multiple_regions, transcripts=input$transcripts_sel), idcol = "origin")
    colnames(transcripts_DF) <- column_names
    transcripts_DF <- transcripts_DF %>%
      filter(origin %in% genoSelect())
    updateSelectizeInput(session, "transcriptList_columnSel", choices=colnames(transcripts_DF))
    shinyalert("Nice!", paste(nrow(transcripts_DF), "transcripts in table"), type = "success")
    transcripts_DATA <- list(transcripts_DF = transcripts_DF, GTF_DF = GTF_DF)
    
    return(transcripts_DATA)

  })
  
  ## transcript specific data (FLEPseq results, gtf and coordinates)
  transcript_data <- eventReactive(input$Submittranscript,{
    req(input$runSelection, input$transcript_sel)
    # GTF gymnastics to handle introns
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
    column_names <- c('origin', column_names)
    tabixed_list <- global$sample_corr$tabix_file
    names(tabixed_list) <- global$sample_corr$genotype
    tabixed_df = setDT(as.list(tabixed_list))
    tryCatch(
      {
        transcript_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files_single_region, transcript=input$transcript_sel), idcol = "origin")
        shinyalert("Nice!", paste(nrow(transcript_DF), "transcripts in table"), type = "success")
      },
      error = function(cond) {shinyalert("OUCH", paste("Error:", cond), type = "error")}
    )
    
    tryCatch(
      expr = {
        colnames(transcript_DF) <- column_names
        transcript_DF <- transcript_DF %>%
          group_by(origin) %>% 
          mutate(sample_id = row_number()) %>% 
          ungroup() %>%
          mutate(Run=basename(global$datapath))%>%
          mutate(across(retention_introns, as.character))%>%
          mutate(retention_introns = replace_na(retention_introns, "none")) %>%
          filter(origin %in% genoSelect())
        
        n_introns <- unique(transcript_DF$mRNA_intron_num)
        if (n_introns>0) {
          intron_cols <- paste0("intron",1:n_introns)
          setDT(transcript_DF)
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
               eventExpr = {input$runSelection},
               handlerExpr = {
                 # Sets paths variables
                 
                 if (!input$runSelection == "") {

                   global$datapath <- input$runSelection
                   tail_files <- list.files(path=file.path(global$datapath,tail_dir), pattern = paste0(tail_ext, "$"),full.names = T)
                   
                   sample_file <- list.files(path=global$datapath, pattern=sample_table, full.names = T)
                   mapping_files <- list.files(path=file.path(global$datapath,mapping_dir), pattern=paste0(mapping_ext, "$"), full.names=T)
                   
                   if (length(tail_files)>0) {
                     global$sample_corr <- fread(sample_file, col.names = c("sample", "genotype"), header = F) %>%
                       mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(sample, index_ext)),
                              tail_file=file.path(global$datapath, tail_dir, paste0(sample, tail_ext)),
                              map_file=file.path(global$datapath, mapping_dir, paste0(sample, mapping_ext)),
                              gene_list=file.path(global$datapath, tail_dir, paste0(sample, tabix_l_ext)))
                     
                     # Check if files are tabixed, and tabix them if not
                     apply(global$sample_corr[,c('tabix_file','tail_file')], 1, function(y) check_if_tabixed(y['tabix_file'],y['tail_file']))
                     
                   }  else {
                     global$sample_corr <- fread(sample_file, col.names = c("sample", "genotype"), header = F) %>%
                       mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(sample, index_ext)),
                              map_file=file.path(global$datapath, mapping_dir, paste0(sample, mapping_ext)),
                              gene_list=file.path(global$datapath, tail_dir, paste0(sample, tabix_l_ext)))
                                          
                   }
                   tabix_files <- paste(basename(global$sample_corr$tabix_file), collapse='\n')
                   shinyalert("Nice!", paste("Successfully added", tabix_files, sep="\n"), type = "success")
                   genes_list=unique(rbindlist(lapply(global$sample_corr$gene_list, fread, header=F)))
                   updateSelectizeInput(session, 'transcript_sel', label=NULL, selected="", choices = genes_list$V1, options = list(create = FALSE), server = TRUE)
                   updateSelectizeInput(session, 'transcripts_sel', label=NULL, selected="", choices = genes_list$V1, options = list(create = FALSE), server = TRUE)
                   
                   
                 }
               }) # end observer runSelection
  
  

  ### OUTPUTS ---------------------------------------------------------------

  ## intronic profile selection (radiobuttons) ----
  output$checkbox_introns <- renderUI(radioButtons('retention_introns', 'Select Retention intron', choices =unique(transcript_data()$transcript_DF$retention_introns), inline = T))

  
  output$sample_table = DT::renderDataTable({
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
  
  ## table with user-selected column (1 transcript) ----
  output$FlepTable_single  <- renderDataTable(filtereddata_single())
  
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
  output$FlepTable_list  <- renderDataTable(filtereddata_list())
  
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
    req(transcript_data()$GTF_DF ),
    options = list(pageLength = 50)
  )
  
  output$GTFtable_list <- renderDataTable(
    req(transcripts_data()$GTF_DF ),
    options = list(pageLength = 50)
  )
  
  ## coverage plot ----
  output$coverage <- renderPlot({
    req(MAP_data())
    cov <- MAP_data() %>%
      group_by(origin, rname) %>%
      summarise(mean_cov=mean(coverage),
                mean_mapq =mean(meanmapq))
    ggplot(cov, aes(x=rname, y=mean_cov, fill=origin)) +
      geom_col(position = "dodge") +
      ggtitle("Mean coverage") +
      ggcustom_theme +
      theme(
        legend.position = c(.05, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid", colour ="black")
      )
  })
  
  ## mapping quality plot ----
  output$mapq <- renderPlot({
    req(MAP_data())
    cov <- MAP_data() %>%
      group_by(origin, rname) %>%
      summarise(mean_cov=mean(coverage),
                mean_mapq =mean(meanmapq))
    ggplot(cov, aes(x=rname, y=mean_mapq, fill=origin)) +
      geom_col(position = "dodge") +
      ggtitle("Mean mapping quality")+
      ggcustom_theme +
      theme(legend.position = "none")
  })
  
  ## polyA bulk distribution for single transcript ----
  output$polyA_distr_transcript <- renderPlot({
    req(transcript_data())
    plot_polya_bulk(transcript_data()$transcript_DF)
  })
  
  ## polyA bulk distribution for multiple transcripts ----
  output$polyA_distr_transcripts <- renderPlot({
    req(transcripts_data())
    plot_polya_bulk(transcripts_data()$transcripts_DF)
    
  })
  
  
  ## gene plot ----
  output$plot_gene <- renderPlot({
    req(transcript_data())
    GTF_DF <- transcript_data()$GTF_DF 
    gtf_gene <- GTF_DF%>% filter(feat_type=="gene")
    gtf_subgene <- GTF_DF%>% filter(feat_type=="subgene")
    parent_start <- gtf_gene$start
    parent_stop <- gtf_gene$end
    title <- paste(gtf_gene$transcript, gtf_gene$ROI)
    ggplot() +
      # plot model gene...
      geom_gene_arrow(data=gtf_gene,
                      aes(xmin = start, xmax = end, y = seqnames,forward=orientation), fill = "white") +
      # ...annotated
      geom_subgene_arrow(data = gtf_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = 1, forward=orientation, fill = feature,
                             xsubmin = start, xsubmax = end), color="black") +
      geom_subgene_label(
        data = gtf_subgene,
        aes(y= seqnames, xsubmin = start, xsubmax = end, label = feat_id),
        min.size = 0
      ) +
      ggtitle(title) +
      theme(legend.position="bottom")

  })
  
  
  ## plot reads (intronic profile) ----
  output$plot_reads <- renderPlot({
    req(filteredintron())
    
    validate(need(!is.null(input$retention_introns), "Please select an intron profile"))
    total_data <- build_intronic_profile_for_plot(filteredintron())

    # plot
    ggplot() +
      # transcripts.....
      geom_rect(data=total_data$df_coords_transcript,
                aes(xmin = start, xmax = end, ymin = sample_id, ymax=sample_id + 0.5 ), fill="grey",alpha=.8)+
      # .... annotated
      geom_rect(data = total_data$df_coords_transcript_subgene_tail,
                         aes(xmin = start, xmax = end, ymin = sample_id , ymax=sample_id +0.5,fill = feature), alpha=.8 ) +
      
      # model gene ...
      geom_rect(data=total_data$df_coords_gene, 
                      aes(xmin = start, xmax = end, ymin = 0, ymax= 0.5), color = "black") +
      # ... annotated
      geom_rect(data = total_data$df_coords_subgene,
                        aes(xmin = start, xmax = end, ymin = 0, ymax= 0.5 , fill = feature)) +
       
      facet_wrap(~origin, ncol = 1) +
      theme(strip.background =element_rect(fill="darkgrey"))+
      theme(strip.text = element_text(colour = 'white')) +
      theme(legend.position="bottom") +
      theme(legend.background = element_rect(linewidth=0.5, linetype="solid")) +
      ggtitle(total_data$gene_name) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    
  }, height = 1500)
  
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot_reads_dblclick, {
    brush <- input$plot_reads_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
 
} # /server

# App calling ------------
shinyApp(ui, server)

