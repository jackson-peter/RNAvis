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
    
    selectInput("selectdir", 'Select File', c("Choose one" = "", list.files(datasets_path, recursive = F, full.names = T))),

    br(),
    h5(HTML("2 | Enter an AGI")),
    textInput("AGI",value = "AT1G01010.1",label = NULL ), #textInput
    actionButton(inputId = "SubmitAGI",
                 label = "Get table",
                 style="color: #fffc300; background-color: #e95420; border-color: #c34113;
                                                         >border-radius: 10px; >border-width: 2px"),
    
    h5(HTML("3 | Select usefull columns")),
    selectizeInput('column_sel', NULL,choices=NULL, selected=NULL, multiple=T),
    actionButton(inputId = "update",
                 label = "Update columns",
                 style="color: #fffc300; background-color: #e95420; border-color: #c34113;
                                                         >border-radius: 10px; >border-width: 2px"),
    
    br(),
    h5(HTML("4 | Export Table")),
    checkboxInput("save_selection", "Only save selected columns?", value = FALSE, width = NULL),
    actionButton("save_df","Export Table",
                 style="color: #fffc300; background-color: #e95420; border-color: #c34113;
                                                         >border-radius: 10px; >border-width: 2px")
    
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
                     mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(barcode, ".sorted.csv.gz")),
                            tail_file=file.path(global$datapath, tail_dir, paste0(barcode, tail_ext)),
                            map_file=file.path(global$datapath, mapping_dir, paste0(barcode, mapping_ext)))
        # Check if files are tabixed, and tabix them if not
        apply(global$barcode_corr[,c('tabix_file','tail_file')], 1, function(y) check_if_tabixed(y['tabix_file'],y['tail_file']))
        
      } 
    }
  )
  

  
  # Mapping data
  MAP_data <- eventReactive(input$selectdir,{
    req(global$barcode_corr)
    print(global$barcode_corr)
    map_files <- global$barcode_corr$map_file
    names(map_files) <- global$barcode_corr$genotype
    mapping_q <- rbindlist(lapply(map_files, fread, col.names=mapping_cols), idcol = "origin")
  })
  
  # AGI_data is the data table of the AGI
  AGI_data <- eventReactive(input$SubmitAGI,{
    req(input$selectdir, input$AGI)
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
        #print(head(AGI_DF))
        updateSelectizeInput(session, "column_sel", choices=colnames(AGI_DF))
      },
      error = function(e){ 
        shinyalert("Erf!", "It seems that there's no data in the files chosen for the requested AGI", type = "warning")
      }
    )
    
    return(AGI_DF)

  })
  
  filtereddata <- eventReactive({
    input$update
    AGI_data()
  },  {
    req(AGI_data())
    if (is.null(input$column_sel) || input$column_sel == "") {
      AGI_data()
    } else {
      print(input$column_sel)
      AGI_data()[ ,colnames(AGI_data()) %in% input$column_sel, with=FALSE]
    } 
  })
  
  observeEvent(AGI_data(), {
    updateSelectInput(session, "select", choices=colnames(AGI_data()))
  })

  output$run_name <- renderText({
    #req(AGI_data())
    req(MAP_data())
    title <- basename(global$datapath)
  })
  
  output$AGI_name <- renderText({
    req(AGI_data())
    title <- unique(AGI_data()$mRNA)
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
  
  output$AGI_dt_t <- renderTable({
    dt_render <- filtereddata()
    validate(
      need(nrow(dt_render) > 0, "No Data to show. Please check you provided a valid AGI.")

    )
    dt_render
  })
  

  output$AGI_dt <- renderDataTable(filtereddata(), options = list(scrollX = TRUE))
  
  output$barcode_geno <- renderDataTable(
    global$barcode_corr,
    options = list(scrollX = TRUE)
    )
  
  observeEvent(AGI_data(),updateSelectizeInput(session, "column", choices=names(AGI_data())))
  
  
  ######## SAVE TABLE
  # Fait débuter le browser de fichiers de sauvegarde de df dans le home.
  shinyDirChoose(
    input,
    'save_dir',
    roots = c(home = './data')
  ) #shinyDirChoose
  
  save_global <- reactiveValues(datapath = getwd())
  save_dir <- reactive(input$save_dir)
  
  output$save_dir <- renderText({ 
    save_global$datapath 
  }) #output$dir
  
  # Met a jour le path du dossier choisi.
  observeEvent(ignoreNULL = TRUE, eventExpr = {input$save_dir},
               handlerExpr = {
                 if (!"path" %in% names(save_dir())) return()
                 home <- normalizePath('./data')
                 save_global$datapath <- file.path(home,
                                                   paste(unlist(save_dir()$path[-1]),
                                                         collapse = .Platform$file.sep))
               }) #observeEvent
  
  # Fenêtre de confirmation de sauvegarde du dataframe.
  observeEvent(input$save_df, {
    default_name <- paste0(input$AGI, "tails.tsv")
    showModal(modalDialog(
      tagList(
        textInput("filename", label = "Filename", placeholder = default_name, value = default_name),
        
        shinyDirButton("save_dir", 'Select where to save the dataframe', 'Please select a folder', FALSE),
        verbatimTextOutput("save_dir", placeholder = TRUE)
      ), 
      title="Save the dataframe as .tsv",
      footer = tagList(actionButton("confirmCreate", "Create"),
                       modalButton("Cancel")
      )
    ))
  })
  
  # Enregistre le dataframe en .tsv.
  observeEvent(input$confirmCreate, { 
    req(input$filename)
    removeModal()
    
    show_modal_spinner(spin = "fingerprint", text=paste("exporting table")) # show the modal window
    if (input$save_selection) {
      write_tsv(filtereddata(), file.path(save_global$datapath, input$filename))
    } else {
      write_tsv(AGI_data(), file.path(save_global$datapath, input$filename))
    }
    
    shinyalert("Nice!", paste("Export successfull", file.path(save_global$datapath, input$filename)), type = "success")
    remove_modal_spinner() # remove it when done
  })
  
} #server

shinyApp(ui, server)