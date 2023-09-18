source("global.R")


#######################################################################################################################
#######################################################################################################################
############################################       FUNCTIONS       ####################################################
#######################################################################################################################
#######################################################################################################################

sort_and_tabix <- function(tail_f) {
  bgzip <- file.path(HTSLIB_PATH,"bgzip")
  tabix <- file.path(HTSLIB_PATH,"tabix")
  
  tail_f_sorted <- gsub(extension, ".sorted.csv", tail_f)
  tail_f_gz <- paste0(tail_f_sorted, ".gz")
  sort_cmd=paste0("{ head -n 1 ",tail_f, " && tail -n +2 ", tail_f, " | sort -k",AGI_col," -k",start_col,"n,",end_col,"n; } > ", tail_f_sorted)
  tabix_cmd <- paste(bgzip, tail_f_sorted, "&&",  tabix, tail_f_gz, "-S 1 -s",AGI_col,"-b",start_col,"-e",end_col)
  
  system(sort_cmd)
  system(tabix_cmd)
}

check_if_tabixed <- function(tabix_file,tail_file, index=".tbi") {
  if (!file.exists(tail_file)) {
    shinyalert("FFFFFFF", paste("THERE IS NO", tail_file), type = "error")
  }

  if (!file.exists(tabix_file) | !file.exists(paste0(tabix_file, index))){
    
    show_modal_spinner(spin = "fingerprint", text=paste("indexing", basename(file.path(tabix_file)))) # show the modal window
    sort_and_tabix(tail_file)
    remove_modal_spinner() # remove it when done

    
  } else {
    shinyalert("Nice!", paste("successfully added", tabix_file), type = "success")
  }
}

read_tabixed_files <- function(file, AGI) {
  dt <- fread(cmd = paste(file.path(HTSLIB_PATH, "tabix"), file, AGI, "-h"))
  return(dt)
}
  

#######################################################################################################################
#######################################################################################################################
############################################       UI       ###########################################################
#######################################################################################################################
#######################################################################################################################
ui <- dashboardPage(
  dashboardHeader(title = "RNA visualization",
                  dropdownMenuOutput("notificationMenu")), #dashboardHeader
  
  dashboardSidebar(
    sidebarMenu( 
      br(),
      h4(HTML("&nbsp; 1 | Select a FLEPseq run folder")),
      shinyDirButton("dir", 'Select a folder', 'Please select a run folder  (use the arrows)'),
      

      br(),
      
      h4(HTML("&nbsp; 2 | Enter an AGI to search for")),
      textInput("AGI",value = "AT1G01010.1",label = NULL ), #textInput
      
      # Select variables to display ----
      uiOutput("checkbox"),
      actionButton(inputId = "SubmitAGI",
                   label = "Submit AGI NOW",
                   style="color: #fffc300; background-color: #e95420; border-color: #c34113;
                                                         >border-radius: 10px; >border-width: 2px"),

      
      br(),
     
      
      h4(HTML("&nbsp; 5 | Export Table")),
      
      actionButton("save_df",
                   "Export Table")
    ) #sidebarMenu
  ), #dashboardSidebar
  
  dashboardBody(
    fluidRow(
      tabsetPanel(
        tabPanel("mRNA statistics",
                 tabsetPanel(
                   title = NULL,
                   # The id lets us use input$tabset1 on the server to find the current tab
                   id = "tabset1",
                   
                   tabPanel("DataTable",
                            dataTableOutput("AGI_dt")
                            
                   ),
                   
                 ) #tabsetPanel
        ) #tabPanel
      ) #tabsetPanel
    ),
    fluidRow(infoBoxOutput("tabset1Selected"))
  )
) #ui


#######################################################################################################################
#######################################################################################################################
############################################       SERVER       #######################################################
#######################################################################################################################
#######################################################################################################################


server <- function(input, output, session) {
  
  
  # Fait débuter le browser de fichiers dans le home.
  shinyDirChoose(
    input,
    'dir',
    roots = c(home = './data')
  ) #shinyDirChoose
  
  
  global <- reactiveValues(datapath = getwd())
  dir <- reactive(input$dir)
  output$dir <- renderText({ 
    global$datapath 
  }) #output$dir
  
  
  # Met a jour le path du dossier choisi.
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$dir},
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("./data")
                 global$datapath <- file.path(home,
                                              paste(unlist(dir()$path[-1]),
                                                    collapse = .Platform$file.sep))
                 global$tail_files <- list.files(path=file.path(global$datapath,tail_dir), pattern = paste0(extension, "$"),full.names = T)
                 global$barcode_file <- file.path(global$datapath, sample_table)
               }) #observeEvent
  

  # Actualise les checkboxes avec les noms des barcodes présents dans le dossier.
  observeEvent(
    input$dir, {
      
      #tail_files=list.files(path=file.path(global$datapath,tail_dir), pattern = paste0(extension, "$"),full.names = T)
      tail_files=global$tail_files

      if (length(tail_files)>0) {

        barcode_corr <- fread(global$barcode_file, col.names = c("barcode", "genotype"), header = F) %>%
          mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(barcode, ".sorted.csv.gz")),
                 tail_file=file.path(global$datapath, tail_dir, paste0(barcode, extension)))
        
        # Check if files are tabixed, and tabix them if not
        apply(barcode_corr[,c('tabix_file','tail_file')], 1, function(y) check_if_tabixed(y['tabix_file'],y['tail_file']))
        
      }
    }
  )
  
  # AGI_data is the data table of the AGI
  AGI_data <- eventReactive(input$SubmitAGI,{
    print("building DF")
    barcode_corr <- fread(global$barcode_file, col.names = c("barcode", "genotype"), header = F) %>%
      mutate(tabix_file= file.path(global$datapath, tail_dir, paste0(barcode, ".sorted.csv.gz")),
             tail_file=file.path(global$datapath, tail_dir, paste0(barcode, extension)))
    
    column_names <- names(fread(cmd = paste('head -n 500', barcode_corr$tail_file[1])))
    column_names <- c('origin', column_names)
    tabixed_list <- barcode_corr$tabix_file
    names(tabixed_list) <- barcode_corr$genotype
    
    
    tabixed_df = setDT(as.list(tabixed_list))
    

    AGI_DF <- rbindlist(lapply(tabixed_df, read_tabixed_files, AGI=input$AGI), idcol = "origin")
    colnames(AGI_DF) <- column_names
    #print(head(AGI_DF))
    return(AGI_DF)
  })
  
  output$AGI_dt <- renderDataTable(AGI_data())
  

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
    showModal(modalDialog(
      tagList(
        textInput("filename", label = "Filename", placeholder = paste0(input$AGI, "_tails.tsv")),
        
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
    
    # Requiert un nom de fichier pour créer la sauvegarde.
    req(input$filename)
    
    # Close the confirmation window.
    removeModal()
    
    show_modal_spinner(spin = "fingerprint", text=paste("exporting table")) # show the modal window
    write_tsv(AGI_data(), file.path(save_global$datapath, input$filename))
    shinyalert("Nice!", paste("Export successfull"), type = "success")
    remove_modal_spinner() # remove it when done
  })
  
  
  
} #server



shinyApp(ui, server)