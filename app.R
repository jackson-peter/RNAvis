source("helper_functions.R")
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

# header ------
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
                   )

# sidebar ------
sidebar <-
  dashboardSidebar(
    hr(),
    sidebarMenu(id = "tabs",
                menuItem("Main Dashboard", tabName = 'dashboard', icon = icon('dashboard', selected=TRUE)),
                menuItem("Plots", tabName = "plots", icon=icon("line-chart")),
                menuItem("Tables", tabName = "tables", icon=icon("table"))
                
    ) #/ sidebarmenu
  ) #/ dashboardsidebar

body <-
  dashboardBody(
  tags$head(includeCSS("www/custom.css")),
  # Sidebar Tab Main dashboard -----
  tabItems(
    # sidebar Tab: dashboard
    tabItem(tabName = "dashboard", 
            # Inputs box ---
            fluidRow(
              box(width = 12, status = "primary", solidHeader = TRUE, title="Input files selector",
                  column(width=12, 
                         selectInput("selectdir", NULL, choices = c("Choose one" = "", flep_runs))), # input dir input field
                  column(width=12,
                         selectizeInput("AGI", inputId = 'AGI', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE))),
                  column(width=3,
                         actionButton(inputId = "SubmitAGI", label = "Get table")),
                  column(width=12,
                         textOutput("rundir"),
                         textOutput("agi")))

            ), 
            # --- Inputs box
            # recap of selected inputs ---
            fluidRow(
              box(width = 12, status = "primary", solidHeader = TRUE, title="Input files details", collapsible = T, collapsed = T,
                  column(width=12, 
                         fluidRow(dataTableOutput("barcode_geno")))
              ), 
              
            ), 
            # ---recap of selected inputs
            
            # plots about run
            fluidRow(splitLayout(cellWidths = c("50%", "50%"), plotOutput("coverage"), plotOutput("mapq")))
            
    ), # /tabitem dashboard
    
    # Sidebar Tab Plots -----
    tabItem(tabName = "plots",
            fluidRow(
              tabBox(title = "Plots",
                     width = NULL,
                     tabPanel(h5("READS"),uiOutput("plotreads.ui"))
                     #tabPanel(h5("READS"),uiOutput("plotreads.ui")),
                     ),
            ),
    ),

    # Sidebar Tab Tables -----
    tabItem(tabName="tables",
            fluidRow(
              tabBox(title = "Tables",
                     width=NULL,
                     # Tab with GFF infos
                     tabPanel(h5("Annotation"), dataTableOutput("GFF_table")),
                     # Tab with AGI table
                     tabPanel(h5("AGI infos"), 
                              fluidRow(
                                selectizeInput('column_sel', NULL,choices=NULL, selected=NULL, multiple=T), # columns selector
                                actionButton(inputId = "update",
                                             label = "Update columns"),
                                checkboxInput("save_selection", "Only save selected columns?", value = FALSE, width = NULL),
                                downloadButton("downloadData", "Download"),
                                uiOutput("AGI_dt_t")
                                       
                              ) # /fluidRow
        
                    ) # /tabPanel

              ) # /tabbox
            ) #/ fluidrow
    ) #/ tabitem

  ) #/ tabitems
) #/ dashboardbody 


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
  
  ### OUTPUTS ---
  text_rundir <- reactiveVal('Please select a directory')
  output$rundir <- renderText({text_rundir()}) #output$rundir
  
  text_agi <- reactiveVal('Please select a gene')
  output$agi <- renderText({text_agi()}) #output$rundir
  
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
  
  output$AGI_dt_t <- renderTable({
    dt_render <- filtereddata()
    validate(
      need(nrow(dt_render) > 0, "No Data to show. Please check you provided a valid AGI.")
    )
    dt_render
  })
  
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
    })
  
  
  
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
    df_coords_gene <- coords_df %>% 
      filter(feat_type=="gene") 
    
    df_coords_gene <- df_coords_gene[!duplicated(df_coords_gene$mRNA, df_coords_gene$retention_introns),]
    
    df_coords_subgene <- coords_df %>% filter(feat_type=="subgene") %>%
      mutate(parent_start=unique(df_coords_gene$start),
             parent_stop=unique(df_coords_gene$end))
    
    df_coords_tails <- coords_df %>% filter(feat_type=="tail") %>%
      mutate(parent_start=unique(df_coords_gene$start),
             parent_stop=unique(end))
    
    
    df_coords_transcript <- coords_df %>% filter(feat_type=="transcript")
    df_coords_transcript_subgene <- df_coords_subgene %>%
      filter(start>=read_start,
             end<=read_end)
    
    df_coords_transcript_subgene_tail <- rbind(df_coords_transcript_subgene, df_coords_tails)
    
    df_coords_transcript_subgene_tail$parent_start <- min(df_coords_transcript_subgene_tail$start, na.rm = T)
    df_coords_transcript_subgene_tail$parent_stop <- max(df_coords_transcript_subgene_tail$end, na.rm=T)
    
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
  
  
  ploutCountRead <- reactive({
    req(AGI_data())
    return(nrow(AGI_data()$AGI_DF))
  })
  
  plotHeight <- reactive(20 * ploutCountRead()) 
  
  output$plotreads.ui <- renderUI({
    plotOutput("plot_reads", height = plotHeight())
  })
  
  
  ### --- OUTPUTS
  
  ### OBSERVERS ---
  
  # Observer for selectdir. 
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$selectdir},
               handlerExpr = {
                 # Sets paths variables
                 if (!input$selectdir == "") {
                   text_rundir(paste("Selected directory:",input$selectdir, sep=" ")) # set new value to reactiveVal
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
                     updateSelectizeInput(session, 'AGI', label=NULL, selected="", choices = genes_list$V1, options = list(create = FALSE), server = TRUE)
                     
                   } 
                 }
               }) # end observer selectdir
  
  # Observer for AGI 
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$AGI},
               handlerExpr = {
                 if (!input$selectdir == "") {
                   text_agi(paste("Selected gene:",input$AGI, sep=" "))}
               })
  
  ### --- OBSERVERS
  
  ### EventReactive (reactive DATASETS) ---
  
  # Mapping data
  MAP_data <- eventReactive(input$selectdir,{
    req(global$barcode_corr)
    map_files <- global$barcode_corr$map_file
    names(map_files) <- global$barcode_corr$genotype
    mapping_q <- rbindlist(lapply(map_files, fread, col.names=mapping_cols), idcol = "origin")
  })
  
  
  # AGI specific data (FLEPseq results, gff and coordinates)
  AGI_data <- eventReactive(input$SubmitAGI,{
    
    req(input$selectdir, input$AGI)
    
    # GFF gymnastics to handle introns
    gff_infos <- read_GFF_file(ref_gff, input$AGI)
    mRNA_df <- gff_infos %>% filter(feature=="mRNA")
    exon_df <- gff_infos %>% filter(feature=="exon")
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
    
    # Building AGI dataset
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
        
        # Building coords df
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
  
  # filters AGI_data with users' selection of columns
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
  
  
  
  ### --- EventReactive (reactive DATASETS)


} #server

### App calling       
shinyApp(ui, server)

