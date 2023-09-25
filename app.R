source("global.R")

#######################################################################################################################
#######################################################################################################################
############################################       UI       ###########################################################
#######################################################################################################################
#######################################################################################################################
header= dashboardHeader(title = "RNAvis", dropdownMenuOutput("notificationMenu"))

sidebar <-   dashboardSidebar(
  sidebarMenu(
    br(),
    h5(HTML("1 | Select a FLEPseq folder")),

    selectInput("selectdir", 'Select File', c("Choose one" = "", flep_runs)),

    br(),
    h5(HTML("2 | Enter an AGI")),
    #textInput("AGI",value = "AT5G09810.1",label = NULL ), #textInput
    selectizeInput("AGI", inputId = 'AGI', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE)),
    actionButton(inputId = "SubmitAGI",
                 label = "Get table"),

    h5(HTML("3 | Select usefull columns")),
    selectizeInput('column_sel', NULL,choices=NULL, selected=NULL, multiple=T),
    actionButton(inputId = "update",
                 label = "Update columns"),

    br(),
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
             fluidRow(uiOutput("AGI_dt_t"))
    ), # tabPanel
    
    
    tabPanel("plot_3p_end",
             #h2(textOutput("AGI_name")),
             fluidRow(plotOutput("plot_3p_end")),
             hr(),
    ),
    
    tabPanel("plot_5p_end",
             #h2(textOutput("AGI_name")),
             fluidRow(plotOutput("plot_5p_end")),
             hr(),
    ),
    
  ), #tabsetPanel
) # /dashboardbody


### UI construction

ui <- dashboardPage(skin="green",
  header, 
  sidebar,
  body
)

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
  
  output$run_name <- renderText({
    req(MAP_data())
    title <- basename(global$datapath)
  })
  
  # AGI_data is the data table of the AGI
  AGI_data <- eventReactive(input$SubmitAGI,{
    req(input$selectdir, input$AGI)
    print(input$AGI)
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
        updateSelectizeInput(session, "column_sel", choices=colnames(AGI_DF))
        
      },
      error = function(e){ 
        shinyalert("Erf!", "It seems that there's no data in the files chosen for the requested AGI", type = "warning")
      }
    )
    print(head(AGI_DF))
    return(AGI_DF)
  })
  
  output$AGI_name <- renderText({
    req(AGI_data())
    title <- unique(AGI_data()$mRNA)
  })
  
  filtereddata <- eventReactive({
    input$update
    AGI_data()
  },  {
    req(AGI_data())
    if (is.null(input$column_sel) || input$column_sel == "") {
      AGI_data()
    } else {
      
      AGI_data()[ ,colnames(AGI_data()) %in% input$column_sel, with=FALSE]
    } 
  })
  
  output$AGI_dt_t <- renderTable({
    dt_render <- filtereddata()
    validate(
      need(nrow(dt_render) > 0, "No Data to show. Please check you provided a valid AGI.")
    )
    dt_render
  })

  output$coverage <- renderPlot({
    # req(AGI_data())
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
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black")
      )
  })

  output$mapq <- renderPlot({
    # req(AGI_data())
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

  output$barcode_geno <- renderDataTable(
    global$barcode_corr,
    options = list(scrollX = TRUE)
    )

  ######## SAVE TABLE
  
  output$downloadData <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      current_datetime <- now()
      formatted_datetime <- format(current_datetime, format = "%Y%m%d_%H%M%S")
      print(formatted_datetime)
      
      paste0(input$AGI, "_", formatted_datetime, ".tsv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      if (input$save_selection) {
        write_tsv(filtereddata(), file)
      }else {
        write_tsv(AGI_data(), file)
      }
    }
  )
  
  # Plot position 3'.
  output$plot_3p_end <- renderPlot({
    req(AGI_data())
    mRNA <- unique(AGI_data()$mRNA)
    data_in <- AGI_data() %>% 
      separate(read_core_id, into=c("readname", "Chr","genome_start_coord", "genome_end_coord"), sep=',', convert = T) %>%
      separate(coords_in_read, into=c("gene_start", "gene_end", "polya_start", "polya_end", "add_start", "add_end", "adapt_start", "adapt_end"), sep=':', convert=T)
    
    ggplot(data_in, aes(genome_start_coord+gene_end)) +
      geom_bar(aes(y = (..count..)/sum(..count..), fill=type)) +
      geom_density() +
      geom_gene_arrow( aes(xmin = unique(mRNA_start), xmax = unique(mRNA_end), y = 0)) +
      geom_gene_label( aes(xmin = unique(mRNA_start), xmax = unique(mRNA_end), y = 0, label = mRNA)) +
      facet_wrap(~origin, ncol=1) +

      labs(title=paste("Counts of reads based on 3' ends of AGI", mRNA, ", n=", nrow(AGI_data()))) +
      ylab("Counts of reads") +
      xlab("3' ends coordinates") + 
      theme(legend.position="bottom",legend.background = element_rect(size=0.5, linetype="solid"))
      
  })
  
  # Plot position 5'.
  output$plot_5p_end <- renderPlot({
    req(AGI_data())
    mRNA <- unique(AGI_data()$mRNA)
    data_in <- AGI_data() %>% 
      separate(read_core_id, into=c("readname", "Chr","genome_start_coord", "genome_end_coord"), sep=',', convert = T) %>%
      separate(coords_in_read, into=c("gene_start", "gene_end", "polya_start", "polya_end", "add_start", "add_end", "adapt_start", "adapt_end"), sep=':', convert=T)
    
    ggplot(data_in, aes(genome_start_coord+gene_start)) +
      geom_bar(aes(y = (..count..)/sum(..count..), fill=type)) +
      geom_density() +
      geom_gene_arrow( aes(xmin = unique(mRNA_start), xmax = unique(mRNA_end), y = 0)) +
      geom_gene_label( aes(xmin = unique(mRNA_start), xmax = unique(mRNA_end), y = 0, label = mRNA)) +
      facet_wrap(~origin, ncol = 1) +
      
      labs(title=paste("Counts of reads based on 5' ends of AGI", mRNA, ", n=", nrow(AGI_data()))) +
      ylab("Counts of reads") +
      xlab("5' ends coordinates") +
      theme(legend.position="bottom",legend.background = element_rect(size=0.5, linetype="solid"))
  
  })

} #server

#######################################################################################################################
#######################################################################################################################
#########################################       App calling       #####################################################
#######################################################################################################################
#######################################################################################################################
shinyApp(ui, server)