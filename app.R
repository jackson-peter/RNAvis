source("helper_functions.R")
source("global.R")

theme_set(theme_bw())

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
                menuItem("Transcripts List", tabName = "transcriptsListTab", icon=icon("list-ul"))
                #menuItem("ReadMe", tabName = "readme", icon=icon("mortar-board")),
                #menuItem("About", tabName = "about", icon = icon("question"))
                
    ), #/ sidebarmenu
    hr()
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
                selectizeInput("selectdir", inputId = 'selectdir', label=NULL, choices = c("Choose one" = "", flep_runs), multiple = FALSE),
                fluidRow(column(dataTableOutput("barcode_geno"), width=11)),
                fluidRow(column(textOutput("rundir"),width=11)),
                fluidRow(column(textOutput("genotypes"),width=11))
                ),# /box
            # plots about run
            box(width = 12, status = "primary", solidHeader = TRUE, title="Run mapping & Coverage", collapsible = T, collapsed=T,
                fluidRow(column(splitLayout(plotOutput("coverage"),
                                            plotOutput("mapq"), cellWidths = c("50%", "50%"))
                                ,width=11))
                ) #/ box

    ), # /tabitem dashboard
    
    # Sidebar tab transcriptSpecific -----
    tabItem(tabName = "transcriptSpecificTab",
            box(width = 12, 
                status = "primary", 
                solidHeader = TRUE, 
                title="Select your Transcript of interest",
                column(width=12,
                       selectizeInput("transcript", inputId = 'transcript', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE))),
                column(width=3,
                       actionButton(inputId = "Submittranscript", label = "Get Data")),
                column(width=12,
                       textOutput("transcript"))
            ), # /box
            tabsetPanel(
              tabPanel(h5("Transcript Overview"),
                       plotOutput("plot_gene"),
                       dataTableOutput("GFF_table")
                       ), # /tabpanel
              tabPanel(h5("Intronic Profiles"),
                       fluidRow(box(width=12, status="primary", solidHeader = T, title="Intronic profile selection", collapsible=T, collapsed=F,
                                    uiOutput("checkbox_introns")
                                    ) # /box
                       ),
                       
                       plotOutput("plot_reads", width = "100%")
                       ), # /tabpanel
              tabPanel(h5("FLEPseq Results"),
                       box(width = 12, status = "primary", solidHeader = TRUE, title="Select a subsets of columns",
                           selectizeInput('column_sel', NULL,choices=NULL, selected=NULL, multiple=T), # columns selector
                           actionButton(inputId = "update_cols",label = "Update columns")),
                       dataTableOutput("transcript_data_col_sel")

                       ), # /tabPanel
              tabPanel(h5("PolyA Distribution"),
                       plotOutput("polyA_distr_transcript")
                       ) # /tabpanel
              ) # / tabsetpanel

            
    ), #/ tabitem transcriptspecific
    
    # Sidebar tab transcriptlist -----
    tabItem(tabName = "transcriptsListTab",
            box(width = 12,
                status = "primary", 
                solidHeader = TRUE, title="Select your Genes of interest",
                column(width=12,
                       textAreaInput("transcripts", "Insert you transcripts ids list here (one by line)", width = "400px")),
                column(width=3,
                       actionButton(inputId = "Submittranscripts", label = "Get Data")),
                column(width=12,
                       textOutput("transcripts"))
            ), # /box
            tabsetPanel(
              tabPanel(h5("transcripts FLEPseq results"),
                       dataTableOutput("multiple_transcripts_table")
                       ), # / tabpanel
              tabPanel(h5("PolyA Distribution"),
                       plotOutput("polyA_distr_transcripts")
                       ) # /tabpanel
                            

            ) # / tabsetpanel

    ) #/ tabitem transcriptlist

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

  ## Mapping data (coverage etc...)
  MAP_data <- eventReactive(input$selectdir,{
    req(global$barcode_corr)
    map_files <- global$barcode_corr$map_file
    names(map_files) <- global$barcode_corr$genotype
    mapping_q <- rbindlist(lapply(map_files, fread, col.names=mapping_cols), idcol = "origin") %>%
      filter(origin %in% genoSelect())
  })
  
  ## Selected Genotypes
  genoSelect <- reactive({
    rows=names(input)[grepl(pattern = "srows_",names(input))]
    paste(unlist(lapply(rows,function(i){
      if(input[[i]]==T){
        index=as.numeric(substr(i,gregexpr(pattern = "_",i)[[1]]+1,nchar(i)))
        
        return(global$barcode_corr$genotype[index])
      }
    })))
    
  })
  
  ## Dataset for list of transcripts
  transcripts_data <- eventReactive(input$Submittranscripts, {
    req(input$selectdir, input$transcripts)
    tabixed_list <- global$barcode_corr$tabix_file
    names(tabixed_list) <- global$barcode_corr$genotype
    tabixed_df = setDT(as.list(tabixed_list))
    column_names <- names(fread(cmd = paste('head -n 1', global$barcode_corr$tail_file[1])))
    column_names <- c('origin', column_names)
    transcripts_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files_multiple_regions, transcripts=input$transcripts), idcol = "origin")
    colnames(transcripts_DF) <- column_names
    transcripts_DF <- transcripts_DF %>%
      filter(origin %in% genoSelect()) 

    return(transcripts_DF)
  })
  
  ## transcript specific data (FLEPseq results, gff and coordinates)
  transcript_data <- eventReactive(input$Submittranscript,{
    req(input$selectdir, input$transcript)
    # GFF gymnastics to handle introns
    gff_infos <- read_GFF_file(ref_gff, input$transcript)
    mRNA_df <- gff_infos %>% filter(feature=="mRNA")
    exon_df <- gff_infos %>% filter(feature=="exon")
    mRNA_GR <- GRanges(mRNA_df$ROI)
    exon_GR <- GRanges(exon_df$ROI)
    # if several exons(and therefore introns)
    if (nrow(exon_df)>1) {
      introns_regions <- as.data.frame(GenomicRanges::setdiff(mRNA_GR, exon_GR, ignore.strand=TRUE)) %>%
        mutate(feature="intron",
               feat_type="subgene",
               ROI=paste0(seqnames, ":", start,"-", end),
               orientation=unique(gff_infos$orientation),
               strand=unique(gff_infos$strand),
               source=unique(gff_infos$source))
      
      GFF_DF <- rbind(gff_infos, introns_regions, fill=T) %>%
        group_by(feature) %>%
        mutate(transcript=input$transcript,
               feat_id=case_when(orientation==0 ~paste0(feature, rev(seq_along(feature))),
                                 orientation==1 ~paste0(feature, seq_along(feature))))%>%
        arrange(start)
    # if only one exon
    } else {
      GFF_DF <- gff_infos %>%
        mutate(transcript=input$transcript,
               feat_id=paste0(feature,"1")) %>%
        arrange(start)
    }
    
    # Building transcript dataset
    column_names <- names(fread(cmd = paste('head -n 1', global$barcode_corr$tail_file[1])))
    column_names <- c('origin', column_names)
    tabixed_list <- global$barcode_corr$tabix_file
    names(tabixed_list) <- global$barcode_corr$genotype
    tabixed_df = setDT(as.list(tabixed_list))
    tryCatch(
      {
        transcript_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files_single_region, transcript=input$transcript), idcol = "origin")
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
          coords_df <- build_coords_df(transcript_DF, GFF_DF, intron_cols) 
        } else {
          coords_df <- build_coords_df(transcript_DF, GFF_DF, vector()) 
        }

        updateSelectizeInput(session, "column_sel", choices=colnames(transcript_DF))
      },
      error = function(e){
        shinyalert("Erf!", paste("Something went wrong", e), type = "error")
        print(e)}
     
    )

    transcript_DATA <- list(transcript_DF = transcript_DF, GFF_DF = GFF_DF, COORDS_DF = coords_df)
    
    return(transcript_DATA)
  })

  ## Dataset filtered by intronic profile
  filteredintron <- reactive({ 
    
    filtered_df<- transcript_data()$COORDS_DF %>% 
        dplyr::filter(retention_introns %in% as.vector(input$retention_introns)) 
    
    return(filtered_df)


  })
  
  ## Dataset column subsets
  filtereddata <- eventReactive(input$update_cols,{
    
    if (is.null(input$column_sel)) {
      filtereddata <- transcript_data()$transcript_DF
    } else {
      filtereddata <- transcript_data()$transcript_DF%>% 
        ungroup()%>%
        select(input$column_sel)
    }
    return(filtereddata)
    
  })
  
  ### OBSERVERS ---------------------------------------------------------------

  ## Observer for selectdir. 
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$selectdir},
               handlerExpr = {
                 # Sets paths variables
                 
                 if (!input$selectdir == "") {
                   text_rundir(paste("Selected directory:",input$selectdir, sep=" ")) # set new value to reactiveVal
                   text_geno(unlist(genoSelect())) # set new value to reactiveVal
                   global$datapath <- input$selectdir
                   
                   tail_files <- list.files(path=file.path(global$datapath,tail_dir), pattern = paste0(tail_ext, "$"),full.names = T)
                   barcode_file <- list.files(path=global$datapath, pattern=sample_table, full.names = T)
                   mapping_files <- list.files(path=file.path(global$datapath,mapping_dir), pattern=paste0(mapping_ext, "$"), full.names=T)
                   
                   if (length(tail_files)>0) {
                     global$barcode_corr <- fread(barcode_file, col.names = c("barcode", "genotype"), header = F) %>%
                       mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(barcode, index_ext)),
                              tail_file=file.path(global$datapath, tail_dir, paste0(barcode, tail_ext)),
                              map_file=file.path(global$datapath, mapping_dir, paste0(barcode, mapping_ext)),
                              gene_list=file.path(global$datapath, tail_dir, paste0(barcode, tabix_l_ext)))
                     
                     # Check if files are tabixed, and tabix them if not
                     apply(global$barcode_corr[,c('tabix_file','tail_file')], 1, function(y) check_if_tabixed(y['tabix_file'],y['tail_file']))
                     tabix_files <- paste(basename(global$barcode_corr$tail_file), collapse='\n')
                     shinyalert("Nice!", paste("Successfully added", tabix_files, sep="\n"), type = "success")
                     genes_list=unique(rbindlist(lapply(global$barcode_corr$gene_list, fread, header=F)))
                     updateSelectizeInput(session, 'transcript', label=NULL, selected="", choices = genes_list$V1, options = list(create = FALSE), server = TRUE)
                  
                     
                   } 
                 }
               }) # end observer selectdir
  
  
  ## Observer for transcript 
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$transcript},
               handlerExpr = {
                 if (!input$selectdir == "") {
                   text_transcript(paste("Selected gene:",input$transcript, sep=" "))
                   
                }
               })
  
  ## Observer for transcripts list 
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$transcripts},
               handlerExpr = {
                 if (!input$selectdir == "") {
                   text_transcripts(paste("Selected genes:",input$transcripts, sep=" "))
                   
                 }
               })



  ### OUTPUTS ---------------------------------------------------------------
  
  ## intronic profile selection (radiobuttons) ----
  output$checkbox_introns <- renderUI(radioButtons('retention_introns', 'Select Retention intron', choices =unique(transcript_data()$transcript_DF$retention_introns), inline = T))

  ## selected directory path ----
  text_rundir <- reactiveVal('Please select a directory')
  output$rundir <- renderText({text_rundir()}) #output$rundir
  
  ## selected transcript (1) ----
  text_transcript <- reactiveVal('Please select a gene')
  output$transcript <- renderText({text_transcript()}) #output$transcript
  
  ##  selected transcripts (list) ----
  text_transcripts <- reactiveVal('Please select genes')
  output$transcripts <- renderText({text_transcripts()}) #output$transcripts
  
  ## selected genotypes ----
  text_geno <- reactiveVal('Please select genotypes')
  output$genotypes <- renderText({text_geno()}) #output$geno
  
  
  ## table of samples ----
  selected_samples <- reactive({input$selected_graph})
  output$barcode_geno = DT::renderDataTable({
    req(global$barcode_corr$genotype)
    #Display table with checkbox buttons
    DT::datatable(cbind(Pick=shinyInput(checkboxInput,"srows_",length(global$barcode_corr$genotype),value=TRUE,width=1), global$barcode_corr),
                  options = list(orderClasses = TRUE,
                                 lengthMenu = c(5, 25, 50),
                                 pageLength = 25 ,
                                 drawCallback= JS(
                                   'function(settings) {
                                     Shiny.bindAll(this.api().table().node());}')
                  ),selection='none',escape=F)
  })
  
  ## table with user-selected column ----
  output$transcript_data_col_sel  <- renderDataTable(filtereddata())
  
  ## GFF table ----
  output$GFF_table <- renderDataTable(
    req(transcript_data()$GFF_DF ),
    options = list(pageLength = 50)
  )
  
  ##  table for multiple transcripts ----
  output$multiple_transcripts_table <- renderDataTable(
    req(transcripts_data()),
    options = list(pageLength = 100)
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
      theme(legend.position = "none")
  })
  
  ## polyA bulk distribution for single transcript ----
  output$polyA_distr_transcript <- renderPlot({
    req(transcript_data())
    ggplot(transcript_data()$transcript_DF, aes(x=polya_length, color=origin)) +
      geom_density()
  })
  
  ## polyA bulk distribution for multiple transcripts ----
  output$polyA_distr_transcripts <- renderPlot({
    req(transcripts_data())
    ggplot(transcripts_data(), aes(x=polya_length, color=origin)) +
      geom_density()
  })
  
  ## poly a "intergenic" distribution for multiple transcripts ----
  output$polyA_distr_transcripts_permrna <- renderPlot({
    req(transcripts_data())
    ggplot(transcripts_data() %>% 
             group_by(mRNA) %>%
             summarise(), aes(x=polya_length, color=origin)) +
      geom_density()
  })
  
  ## gene plot ----
  output$plot_gene <- renderPlot({
    req(transcript_data())
    GFF_DF <- transcript_data()$GFF_DF 
    gff_gene <- GFF_DF%>% filter(feat_type=="gene")
    gff_subgene <- GFF_DF%>% filter(feat_type=="subgene")
    parent_start <- gff_gene$start
    parent_stop <- gff_gene$end
    title <- paste(gff_gene$transcript, gff_gene$ROI)
    ggplot() +
      # plot model gene...
      geom_gene_arrow(data=gff_gene,
                      aes(xmin = start, xmax = end, y = seqnames,forward=orientation), fill = "white") +
      # ...annotated
      geom_subgene_arrow(data = gff_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = 1, forward=orientation, fill = feature,
                             xsubmin = start, xsubmax = end), color="black") +
      geom_subgene_label(
        data = gff_subgene,
        aes(y= seqnames, xsubmin = start, xsubmax = end, label = feat_id),
        min.size = 0
      ) +
      ggtitle(title) +
      theme(legend.position="bottom") +
      theme_void()

  })
  
  ## plot for transcript by  selected intronic profile ----
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
                      aes(xmin = start, xmax = end, ymin = 0, ymax= 0.5), fill = "black") +
      # ... annotated
      geom_rect(data = total_data$df_coords_subgene,
                        aes(xmin = start, xmax = end, ymin = 0, ymax= 0.5 , fill = feature)) +
       
      facet_wrap(~origin, ncol = 1) +
      theme_classic() +
      theme(legend.position="bottom") +
      theme(legend.background = element_rect(linewidth=0.5, linetype="solid")) +
      ggtitle(total_data$gene_name) 
    
  }, height=1500)
 
} # /server

# App calling ------------
shinyApp(ui, server)

