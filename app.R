source("global.R")
source("ui_helper_funct.R")
theme_set(theme_bw())

###
# 1
###
#######################################################################################################################
#######################################################################################################################
############################################       UI       ###########################################################
#######################################################################################################################
#######################################################################################################################
header <- 
  dashboardHeader( title = HTML("RNAvis"),
                   disable = F,
                   dropdownMenuCustom( type = 'message',
                                       customSentence = "ebebebebedb",
                                       messageItem(
                                         from = "jackson.peter",#'Feedback and suggestions',
                                         message =  "",#paste0("jackson.peter@ibmp-cnrs.unistra.fr" ),
                                         icon = icon("envelope"),
                                         href = "mailto:jackson.peter@ibmp-cnrs.unistra.fr"),
                                       icon = icon('comment'))
                   )


sidebar <-
  dashboardSidebar(
    sidebarMenu(
      id="sidebar",
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

    tabPanel("GFF_gene",
             fluidRow(dataTableOutput("GFF_table")),
             
    ),
    tabPanel("Reads viewer",
             fluidRow(plotOutput("plot_reads", height=8000)),
            
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

        if (n_introns>0) {
          intron_cols <- paste0("intron",1:n_introns)
          setDT(AGI_DF)
          AGI_DF <- AGI_DF[, paste(intron_cols) := lapply(paste(intron_cols), getIntronsRetained, retained_introns = as.character(retention_introns)), by = retention_introns]
          coords_df <- build_coords_df(AGI_DF, GFF_DF, intron_cols) 

        } else {
          coords_df <- build_coords_df(AGI_DF, GFF_DF, NA) 
        }

        coords_df <- coords_df %>%
          group_by( chr, read_start, read_end, retention_introns, polya_length, additional_tail) %>%
          arrange(read_end-read_start,retention_introns, polya_length, str_length(additional_tail)) %>%
          mutate(ID = cur_group_id())
        
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
  

  output$plot_reads <- renderPlot({
    req(AGI_data())
    coords_df <- AGI_data()$COORDS_DF 
    print(coords_df)
    df_coords_gene <- coords_df %>% 
      filter(feat_type=="gene") 
    
    df_coords_gene <- df_coords_gene[!duplicated(df_coords_gene$mRNA, df_coords_gene$retention_introns),]
    
    df_coords_subgene <- coords_df %>% filter(feat_type=="subgene") %>%
      mutate(parent_start=unique(df_coords_gene$start),
             parent_stop=unique(df_coords_gene$end))
    
    df_coords_tails <- coords_df %>% filter(feat_type=="tail") 

    df_coords_tails <- df_coords_tails%>%
      mutate(parent_start=unique(df_coords_gene$start),
             parent_stop=unique(end))
    
    df_coords_transcript <- coords_df %>% filter(feat_type=="transcript")

    df_coords_transcript_subgene <- df_coords_subgene %>%
      filter(start>=read_start,
             end<=read_end)
    df_coords_transcript_subgene_tail <- rbind(df_coords_transcript_subgene, df_coords_tails)

    df_coords_transcript_subgene_tail$parent_start <- min(df_coords_transcript_subgene_tail$start)
    df_coords_transcript_subgene_tail$parent_stop <- max(df_coords_transcript_subgene_tail$end)
    
    
    
    ggplot() +
      # plot all transcripts
      geom_gene_arrow(data=df_coords_transcript,
                      aes(xmin = start, xmax = end, y = ID), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
      # and their sub elements
      geom_subgene_arrow(data = df_coords_transcript_subgene_tail,
                         #aes(xmin = parent_start, xmax = parent_stop, y = ID, fill = feature,
                         aes(xmin = parent_start, xmax = parent_stop, y = ID, fill = feature,
                             xsubmin = start, xsubmax = end), color="black", alpha=.4, arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +

      # plot model gene...
      geom_gene_arrow(data=df_coords_gene, 
                      aes(xmin = start, xmax = end, y = -2), fill = "white", arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
      # ...annotated
      geom_subgene_arrow(data = df_coords_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = -2, fill = feature,
                             xsubmin = start, xsubmax = end), color="black", arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +

      facet_grid(origin~retention_introns) 
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