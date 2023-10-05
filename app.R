source("global.R")
theme_set(theme_bw())

###
# 1
###
#######################################################################################################################
#######################################################################################################################
############################################       UI       ###########################################################
#######################################################################################################################
#######################################################################################################################
header= dashboardHeader(title = "RNAvis", dropdownMenuOutput("notificationMenu"))

sidebar <-   dashboardSidebar(
  sidebarMenu(
    h5(HTML("1 | Select a FLEPseq folder")),

    selectInput("selectdir", 'Select a directory', c("Choose one" = "", flep_runs)),


    h5(HTML("2 | Enter an AGI")),
    selectizeInput("AGI", inputId = 'AGI', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE)),
    actionButton(inputId = "SubmitAGI",
                 label = "Get table"),

    h5(HTML("3 | Select usefull columns")),
    selectizeInput('column_sel', NULL,choices=NULL, selected=NULL, multiple=T),
    actionButton(inputId = "update",
                 label = "Update columns"),


    h5(HTML("4 | Export Table")),
    checkboxInput("save_selection", "Only save selected columns?", value = FALSE, width = NULL),
    downloadButton("downloadData", "Download")

  ) #sidebarMenu
)

body <- dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  
  tabsetPanel(
    tabPanel("FLEPseq Runs",
             h2('FLEPseq Runs infos'),
             fluidRow(
               dataTableOutput("flepRunTable")
             ),

             
    ), #tabPanel
    tabPanel("General mRNA statistics",
             h2(textOutput("run_name")),
             fluidRow(
               dataTableOutput("barcode_geno"),
               hr(),
             ),
             fluidRow(
               splitLayout(cellWidths = c("50%", "50%"), plotOutput("coverage"), plotOutput("mapq")),
               hr(),
             ), 

    ), #tabPanel

    tabPanel("Specific AGI",
             h2(textOutput("AGI_name")),
             #fluidRow(uiOutput("AGI_dt_t"))
             fluidRow(uiOutput("AGI_dt_t"))
    ), # tabPanel
    
    tabPanel("plot_5p_end",
             fluidRow(plotOutput("plot_gene5p"),
                      plotOutput("plot_5p_end")),

    ),
    tabPanel("plot_3p_end",
             fluidRow(plotOutput("plot_gene3p"),
                      plotOutput("plot_3p_end")),

    ),
    tabPanel("GFF_gene",
             fluidRow(dataTableOutput("GFF_table")),
             
    ),
    tabPanel("reads",
             fluidRow(plotOutput("plot_reads", height=1600)),
            
    )
    
  ), #tabsetPanel
) # /dashboardbody

### UI construction
ui <- dashboardPage(skin="green",
  header, 
  sidebar,
  body
)

###
# 2
###
#######################################################################################################################
#######################################################################################################################
############################################       SERVER       #######################################################
#######################################################################################################################
#######################################################################################################################
server <- function(input, output, session) {
  
  global <- reactiveValues(datapath = getwd())
  output$dir <- renderText({global$datapath}) #output$dir
    
  # Met a jour le path du dossier choisi, associe les barcodes aux génotypes
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$selectdir},
               handlerExpr = {
                 global$datapath <- input$selectdir
                 
                 tail_files <- list.files(path=file.path(global$datapath,tail_dir), pattern = paste0(tail_ext, "$"),full.names = T)
                 barcode_file <- file.path(global$datapath, sample_table)
                 mapping_files <- list.files(path=file.path(global$datapath,mapping_dir), pattern=paste0(mapping_ext, "$"), full.names=T)
                 
                 
                 if (length(tail_files)>0) {
                   global$barcode_corr <- fread(barcode_file, col.names = c("barcode", "genotype"), header = F) %>%
                     mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(barcode, index_ext)),
                            tail_file=file.path(global$datapath, tail_dir, paste0(barcode, tail_ext)),
                            map_file=file.path(global$datapath, mapping_dir, paste0(barcode, mapping_ext)),
                            gene_list=file.path(global$datapath, tail_dir, paste0(barcode, tabix_l_ext)))
        # Check if files are tabixed, and tabix them if not
        apply(global$barcode_corr[,c('tabix_file','tail_file')], 1, function(y) check_if_tabixed(y['tabix_file'],y['tail_file']))
        
        genes_list=unique(rbindlist(lapply(global$barcode_corr$gene_list, fread, header=F)))
        #updateSelectizeInput(session, "AGI_list", choices=genes_list)
        updateSelectizeInput(session, 'AGI', label=NULL, selected="", choices = genes_list$V1, options = list(create = FALSE), server = TRUE)
        
      } 
    }
  )
  
  # Mapping data
  MAP_data <- eventReactive(input$selectdir,{
    req(global$barcode_corr)
    map_files <- global$barcode_corr$map_file
    names(map_files) <- global$barcode_corr$genotype
    mapping_q <- rbindlist(lapply(map_files, fread, col.names=mapping_cols), idcol = "origin")
    
  })
  
  output$flepRunTable <- renderDataTable({runs_infos})
  
  output$run_name <- renderText({
    req(MAP_data())
    title <- basename(global$datapath)
  })
  
  AGI_data <- eventReactive(input$SubmitAGI,{
    req(input$selectdir, input$AGI)
    gff_infos <- read_GFF_file(ref_gff, input$AGI)
    
    mRNA_df <- gff_infos%>%
      filter(feature=="mRNA")
    exon_df <- gff_infos%>%
      filter(feature=="exon")
    mRNA_GR <- GRanges(mRNA_df$ROI)
    exon_GR <- GRanges(exon_df$ROI)
    if (nrow(exon_df)>1) {
      introns_regions <- as.data.frame(setdiff(mRNA_GR, exon_GR, ignore.strand=TRUE)) %>%
        mutate(feature="intron",
               feat_type="subgene",
               ROI=paste0(seqnames, ":", start,"-", end),
               orientation=unique(gff_infos$orientation),
               strand=unique(gff_infos$strand),
               source=unique(gff_infos$source))
      
      GFF_DF <- rbind(gff_infos, introns_regions, fill=T) %>%
        group_by(feature) %>%
        mutate(AGI=input$AGI,
               feat_id=paste0(feature, seq_along(feature)))%>%
        arrange(start)
    } else {
      GFF_DF <- gff_infos %>%
        mutate(AGI=input$AGI) %>%
        arrange(start)
    }
    
    
    column_names <- names(fread(cmd = paste('head -n 1', global$barcode_corr$tail_file[1])))
    column_names <- c('origin', column_names)
    tabixed_list <- global$barcode_corr$tabix_file
    names(tabixed_list) <- global$barcode_corr$genotype
    tabixed_df = setDT(as.list(tabixed_list))
    AGI_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files, AGI=input$AGI), idcol = "origin")
    tryCatch(
      expr = {
        
        colnames(AGI_DF) <- column_names
        AGI_DF <- AGI_DF %>% mutate(Run=basename(global$datapath))
        n_introns <- unique(AGI_DF$mRNA_intron_num)
        draw_sequence <- c()
        if (n_introns>0) {
          intron_cols <- paste0("intron",1:n_introns)
          # REMOVE THOSE LINES?
          exon_lengths = GFF_DF%>%filter(feature=="exon") %>% pull(width)
          intron_lengths = GFF_DF%>%filter(feature=="intron") %>% pull(width)
          subgene_lengths=paste(cumsum(c(exon_lengths,intron_lengths)[order(c(seq_along(exon_lengths),seq_along(intron_lengths)))]), collapse=':')
          # \REMOVE?
          
          
          setDT(AGI_DF)
          print("EEDGEDGDEG")
          print(intron_cols)

          print(AGI_DF)
          AGI_DF <- AGI_DF[, paste(intron_cols) := lapply(paste(intron_cols), getIntronsRetained, retained_introns = as.character(retention_introns)), by = retention_introns]
          print("EEDGEDGDEG222222222")
          print(AGI_DF)

          coords_df <- build_coords_df(AGI_DF, GFF_DF, intron_cols) %>%
            group_by(read_core_id) %>%
            mutate(ID = cur_group_id())
          

        } 
        updateSelectizeInput(session, "column_sel", choices=colnames(AGI_DF))
        
      },
      error = function(e){ 
        print(e)
        stop(e)
        shinyalert("Erf!", "It seems that something went wrong", type = "warning")
      }
    )
    
    AGI_DATA <- list(AGI_DF = AGI_DF, GFF_DF = GFF_DF, COORDS_DF = coords_df)
    
    return(AGI_DATA)
    
  })
  

  output$barcode_geno <- renderDataTable(
    global$barcode_corr,
    options = list(scrollX = TRUE)
  )
  
  output$GFF_table <- renderDataTable(
    req(AGI_data()$GFF_DF)
  )
  
  output$AGI_name <- renderText({
    req(AGI_data())
    title <- unique(AGI_data()$AGI_DF$mRNA)
  })
  
  filtereddata <- eventReactive({
    input$update
    AGI_data()$AGI_DF
  },  {
    req(AGI_data())
    if (is.null(input$column_sel) || input$column_sel == "") {
      AGI_data()$AGI_DF
    } else {
      
      AGI_data()$AGI_DF[ ,colnames(AGI_data()$AGI_DF) %in% input$column_sel, with=FALSE]
    } 
  })
  
  output$AGI_dt_t <- renderTable({
    dt_render <- filtereddata()
    validate(
      need(nrow(dt_render) > 0, "No Data to show. Please check you provided a valid AGI.")
    )
    dt_render
  })

  ######## SAVE TABLE
  output$downloadData <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      current_datetime <- now()
      formatted_datetime <- format(current_datetime, format = "%Y%m%d_%H%M%S")
      paste0(input$AGI, "_", formatted_datetime, ".tsv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      if (input$save_selection) {
        write_tsv(filtereddata(), file)
      }else {
        write_tsv(AGI_data()$AGI_DF, file)
      }
    }
  )
  
###
# 2.1
###
### PLOTS
  
  # COVERAGE
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
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="black")
      )
  })
  
  # MAPPING QUALITY
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
  

  # Plot position 3'.
  output$plot_3p_end <- renderPlot({
    req(AGI_data())
    mRNA <- unique(AGI_data()$AGI_DF$mRNA)
    data_in <- AGI_data()$AGI_DF %>% 
      separate(read_core_id, into=c("readname", "Chr","genome_start_coord", "genome_end_coord"), sep=',', convert = T) %>%
      separate(coords_in_read, into=c("gene_start", "gene_end", "polya_start", "polya_end", "add_start", "add_end", "adapt_start", "adapt_end"), sep=':', convert=T)
    
    ggplot(data_in, aes(genome_start_coord+gene_end)) +
      geom_bar(aes(y = (..count..)/sum(..count..), fill=type)) +
      geom_density() +
      facet_wrap(~origin, ncol=1) +
      labs(title=paste("Counts of reads based on 3' ends of AGI", mRNA, ", n=", nrow(AGI_data()$AGI_DF))) +
      ylab("Counts of reads") +
      xlab("3' ends coordinates") +
      theme(legend.position="bottom",
            legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="black"))
      
  })
  
  output$plot_gene5p <- renderPlot({
    req(AGI_data())
    gff_gene <- AGI_data()$GFF_DF %>% filter(feat_type=="gene")
    gff_subgene <- AGI_data()$GFF_DF %>% filter(feat_type=="subgene") %>%
      mutate(parent_start=unique(gff_gene$start),
             parent_stop=unique(gff_gene$end))
    
    ggplot(gff_gene, aes(xmin = start, xmax = end, y = source, forward=orientation)) +
      #facet_wrap(~ molecule, scales = "free", ncol = 1) +
      geom_gene_arrow(fill = "white") +
      geom_subgene_arrow(data = gff_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = source, fill = feature,
                             xsubmin = start, xsubmax = end), color="black", alpha=.7) +
      theme_genes()
  })
  
  output$plot_gene3p <- renderPlot({
    req(AGI_data())
    gff_gene <- AGI_data()$GFF_DF %>% filter(feat_type=="gene")
    gff_subgene <- AGI_data()$GFF_DF %>% filter(feat_type=="subgene") %>%
      mutate(parent_start=unique(gff_gene$start),
             parent_stop=unique(gff_gene$end))
    

    
    ggplot(gff_gene, aes(xmin = start, xmax = end, y = source, forward=orientation)) +
      #facet_wrap(~ molecule, scales = "free", ncol = 1) +
      geom_gene_arrow(fill = "white") +
      geom_subgene_arrow(data = gff_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = source, fill = feature,
                             xsubmin = start, xsubmax = end), color="black", alpha=.7) +
      theme_genes()
  })
  
  output$plot_reads <- renderPlot({
    req(AGI_data())
    coords_df <- AGI_data()$COORDS_DF
    #reads=c("e8b23c2d-e2dd-45bd-983c-cacb15bc9cf7", "a00abf0f-75ab-4bfd-85be-9dac7abbf328", "b308939f-4e77-4be6-9c3b-f493bb3dbb93")
    #df <- AGI_data()$COORDS_DF %>% filter(read_id %in% reads)
    df_coords_gene <- AGI_data()$COORDS_DF %>% 
      filter(feat_type=="gene") 
    
    df_coords_gene <- df_coords_gene[!duplicated(df_coords_gene$mRNA, df_coords_gene$retention_introns),]
    
    
    df_coords_subgene <- AGI_data()$COORDS_DF %>% filter(feat_type=="subgene") %>%
      mutate(parent_start=unique(df_coords_gene$start),
             parent_stop=unique(df_coords_gene$end))
    
    df_coords_transcript <- AGI_data()$COORDS_DF %>% filter(feat_type=="transcript")
    # transcript_start <- unique(df_coords_transcript$start)
    # transcript_end <- unique(df_coords_transcript$end)
    # 
    df_coords_transcript_subgene <- df_coords_subgene %>%
      filter(start>=read_start,
             end<=read_end)
   
    
    
    print("GENE")
    print(df_coords_gene)
    write_tsv(df_coords_gene, file="~/DATA/testGENE.tsv")
    print("SUBGENE")
    print(df_coords_subgene)
    write_tsv(df_coords_subgene, file="~/DATA/testSUBGENE.tsv")
    print(df_coords_transcript)
    write_tsv(df_coords_transcript, file="~/DATA/testTRANSCRIPT.tsv")
    

    ggplot() +
      geom_gene_arrow(data=df_coords_transcript,
                      aes(xmin = start, xmax = end, y = ID), fill = "white") +
      
      geom_subgene_arrow(data = df_coords_transcript_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = ID, fill = feature,
                             xsubmin = start, xsubmax = end), color="black", alpha=.4) +
      
      
      geom_gene_arrow(data=df_coords_gene, 
                      aes(xmin = start, xmax = end, y = -1), fill = "white") +
      geom_subgene_arrow(data = df_coords_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = -1, fill = feature,
                             xsubmin = start, xsubmax = end), color="black") +
      
      
      
      
      facet_wrap(~retention_introns, ncol = 1) 
  })
  
  # Plot position 5'.
  output$plot_5p_end <- renderPlot({
    req(AGI_data())
    mRNA <- unique(AGI_data()$AGI_DF$mRNA)
    data_in <- AGI_data()$AGI_DF %>% 
      separate(read_core_id, into=c("readname", "Chr","genome_start_coord", "genome_end_coord"), sep=',', convert = T) %>%
      separate(coords_in_read, into=c("gene_start", "gene_end", "polya_start", "polya_end", "add_start", "add_end", "adapt_start", "adapt_end"), sep=':', convert=T)
    
    
    ggplot(data_in, aes(genome_start_coord+gene_start)) +
      geom_bar(aes(y = (..count..)/sum(..count..), fill=type)) +
      geom_density()+
      facet_wrap(~origin, ncol = 1) +
      labs(title=paste("Counts of reads based on 5' ends of AGI", mRNA, ", n=", nrow(AGI_data()$AGI_DF))) +
      ylab("Counts of reads") +
      xlab("5' ends coordinates") +
      theme(legend.position="bottom",
            legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="black"))
  
  })
  
} #server

###
# 3
###
#######################################################################################################################
#######################################################################################################################
#########################################       App calling       #####################################################
#######################################################################################################################
#######################################################################################################################
shinyApp(ui, server)