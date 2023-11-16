source("helper_functions.R")
source("global.R")

theme_set(theme_bw())

############################################       UI       ###########################################################
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
                   ) #/ dashboardheader

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


# body ------
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
                         fluidRow(DTOutput("barcode_geno"))),
                  column(width=12,
                         selectizeInput("AGI", inputId = 'AGI', label = NULL, choices = NULL, selected = NULL, multiple = FALSE, options = list(create = FALSE))),
                  column(width=3,
                         actionButton(inputId = "SubmitAGI", label = "Validate Parameters")),
                  column(width=12,
                         textOutput("rundir"),
                         textOutput("genotypes"),
                         textOutput("agi"),
                         textOutput("nb_reads")))

            )
            
            
    ), # /tabitem dashboard
    
    # Sidebar Tab Plots -----
    tabItem(tabName = "plots",
            fluidRow(tabBox(title = "Plots",
                            width = NULL,
                            tabPanel(h5("Run Infos"),
                                     # plots about run
                                     fluidRow(
                                       box(width = 12, status = "primary", solidHeader = TRUE, title="Run mapping & Coverage", collapsible = T, collapsed=F,
                                                  column(width=12, 
                                                         splitLayout(cellWidths = c("50%", "50%"), plotOutput("coverage"), plotOutput("mapq")))
                                           ) #/ box
                                       ) #/ fluidRow
                                     ),
                            tabPanel(h5("Intron Retention"),
                                     fluidRow(
                                       box(width=12, status="primary", solidHeader = T, title="Intronic profile selection", collapsible=T, collapsed=F,
                                         column(width=12,
                                                uiOutput("checkbox_introns"))
                                         ) # /box
                                     ), #/ fluidRow
                                     uiOutput("plotreads.ui")),
                            tabPanel(h5("PolyA Site"),plotOutput("plot_polyA_site"))
                            ) # / tabbox
            ) #/ fluidRow
            
    ),

    # Sidebar Tab Tables -----
    tabItem(tabName="tables",
            fluidRow(
              tabBox(title = "Tables",
                     width=NULL,

                     # Tab with AGI table
                     tabPanel(h5("AGI infos"), 
                              fluidRow(
                                selectizeInput('column_sel', NULL,choices=NULL, selected=NULL, multiple=T), # columns selector
                                actionButton(inputId = "update",
                                             label = "Update columns"),
                                
                                DTOutput("AGI_dt_t", height = "auto", fill = T)
                                #uiOutput("AGI_dt_t")
                                       
                              ) # /fluidRow
        
                    ), # /tabPanel

                    # Tab with GFF infos
                    tabPanel(h5("Annotation"), dataTableOutput("GFF_table"))

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

############################################       SERVER       #######################################################


server <- function(input, output, session) {
  # This function stops the execution and fixes the bug requiring R termination when closing app
  session$onSessionEnded(function() {
    stopApp()
  })
  
  global <- reactiveValues(datapath = getwd())
  
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
               feat_id=case_when(orientation==0 ~paste0(feature, rev(seq_along(feature))),
                                 orientation==1 ~paste0(feature, seq_along(feature))))%>%
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
        AGI_DF <- AGI_DF %>% 
          mutate(Run=basename(global$datapath))%>%
          mutate(across(retention_introns, as.character))%>%
          mutate(retention_introns = replace_na(retention_introns, "none")) %>%
          filter(origin %in% genoSelect())
        
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
          group_by(read_id, origin, chr, read_start, read_end, retention_introns, polya_length, additional_tail) %>%
          arrange(retention_introns, read_end-read_start, polya_length, str_length(additional_tail)) %>%
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
  filteredcols <- eventReactive({input$update},{
    req(AGI_data())
    if (is.null(input$column_sel) || input$column_sel == "") {
      return(AGI_data()$AGI_DF)
    } else {
      return(AGI_data()$AGI_DF[ ,colnames(AGI_data()$AGI_DF) %in% input$column_sel, with=FALSE])
    } 
  })
  
  filteredintron <- reactive({ 
    
    if (is.null(input$retention_introns)) {
      return(NULL)
    }  

    filtered_df<- AGI_data()$COORDS_DF %>% 
      filter(retention_introns %in% as.vector(input$retention_introns)) %>%
      arrange(origin) %>%
      group_by(ID, origin)

    return(filtered_df)


  })
  
  ### --- EventReactive 
  
  ### OBSERVERS ---
  observe({
    text_geno(unlist(genoSelect()))
  })
  
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
                   text_agi(paste("Selected gene:",input$AGI, sep=" "))
                   
                }
               })
  
  
  ### --- OBSERVERS

  
  ### OUTPUTS ---
  

  output$checkbox_introns <- renderUI(radioButtons('retention_introns', 'Select Retention intron', choices =unique(AGI_data()$AGI_DF$retention_introns), inline = T))

  text_rundir <- reactiveVal('Please select a directory')
  output$rundir <- renderText({text_rundir()}) #output$rundir
  
  text_agi <- reactiveVal('Please select a gene')
  output$agi <- renderText({text_agi()}) #output$agi
  
  text_geno <- reactiveVal('Please select genotypes')
  output$genotypes <- renderText({text_geno()}) #output$geno
  
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
  
  observe({
    text_geno(genoSelect() )
  })
  
  genoSelect <- reactive({
    rows=names(input)[grepl(pattern = "srows_",names(input))]


    paste(unlist(lapply(rows,function(i){
      if(input[[i]]==T){
        index=as.numeric(substr(i,gregexpr(pattern = "_",i)[[1]]+1,nchar(i)))
        #print(global$barcode_corr$genotype[index])

        return(global$barcode_corr$genotype[index])
      }
    })))
    
  })

  
  
  output$GFF_table <- renderDataTable(
    req(AGI_data()$GFF_DF)
  )
  
  output$AGI_name <- renderText({
    req(AGI_data())
    title <- unique(AGI_data()$AGI_DF$mRNA)
  })
  

  output$AGI_dt_t <- renderDT({
    
    req(AGI_data())
    if (is.null(input$column_sel) || input$column_sel == "") {
      dt_render <- AGI_data()$AGI_DF
    } else {
      dt_render <- filteredcols()
    }
    
    #saveRDS(object = aa.dt(), file = "debugging/aa_dt.rds")
    nlines <- nrow(dt_render)
    current_datetime <- now()
    formatted_datetime <- format(current_datetime, format = "%Y%m%d_%H%M%S")
    outfile <- paste0(input$AGI, "_", formatted_datetime, ".tsv")
    datatable(dt_render, 
              caption = NULL,
              rownames = FALSE, 
              options = list(dom = "Bfrtip", 
                             pageLength = nlines,
                             #deferRender = TRUE,
                             #scrollY = 400,
                             #scroller = TRUE,
                             buttons = list(list(extend = "copy", header = TRUE, title = NULL), 
                                            list(extend = "csv", header = TRUE, filename = outfile))
                             
              ), 
              
              extensions = c("Buttons", "Scroller"),
              
              fillContainer = T)

    
  }, server = FALSE)
  
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
    
    coords_df <- filteredintron()
    validate(need(!is.null(input$retention_introns), "Please select a data set"))
    
    df_coords_gene <- coords_df %>% 
      filter(feat_type=="gene")
    
    genotypes <- unique(coords_df$origin)
    gene_name <- unique(df_coords_gene$mRNA)
    
    df_coords_gene <- df_coords_gene[!duplicated(df_coords_gene$mRNA, by=c("retention_introns", "origin")),]
    df_coords_gene <- df_coords_gene[rep(seq_len(nrow(df_coords_gene)), each = length(genotypes)), ]
    
    df_coords_gene$origin=genotypes
    
    df_coords_subgene <- coords_df %>% filter(feat_type=="subgene") %>%
      mutate(parent_start=unique(df_coords_gene$start),
             parent_stop=unique(df_coords_gene$end))

    df_coords_tails <- coords_df %>% filter(feat_type=="tail") %>%
      mutate(parent_start=unique(df_coords_gene$start),
             parent_stop=unique(end))
    
    df_coords_transcript <- coords_df %>% 
      filter(feat_type=="transcript") 
    
    df_coords_transcript_subgene <- df_coords_subgene %>%
      filter(start>=read_start,
             end<=read_end)
    
    df_coords_transcript_subgene_tail <- rbind(df_coords_transcript_subgene, df_coords_tails)
    
    df_coords_transcript_subgene_tail$parent_start <- min(df_coords_transcript_subgene_tail$start, na.rm = T)
    df_coords_transcript_subgene_tail$parent_stop <- max(df_coords_transcript_subgene_tail$end, na.rm=T)
    
    #show_modal_spinner(spin = "fingerprint", text="Rendering plot") # show the modal window
    limy_df <- df_coords_transcript %>%
      group_by(origin) %>%
      summarise(nb_ID = n_distinct(ID))
    
    ggplot() +
      # plot all transcripts
      geom_gene_arrow(data=df_coords_transcript,
                      aes(xmin = start, xmax = end, y = -1*ID, forward=orientation), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
      # and their sub elements
      geom_subgene_arrow(data = df_coords_transcript_subgene_tail,
                         #aes(xmin = parent_start, xmax = parent_stop, y = ID, fill = feature,
                         aes(xmin = parent_start, xmax = parent_stop, y = -1*ID,forward=orientation, fill = feature,
                             xsubmin = start, xsubmax = end), color="black", alpha=.4, arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
      
      # plot model gene...
      geom_gene_arrow(data=df_coords_gene, 
                      aes(xmin = start, xmax = end, y = 1,forward=orientation), fill = "white", arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
      # ...annotated
      geom_subgene_arrow(data = df_coords_subgene,
                         aes(xmin = parent_start, xmax = parent_stop, y = 1, forward=orientation, fill = feature,
                             xsubmin = start, xsubmax = end), color="black", arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
      
      #facet_grid(origin~retention_introns)
      facet_wrap(~origin, ncol = 1) +
      theme(legend.position="bottom") +
      theme(legend.background = element_rect(size=0.5, linetype="solid")) +
      ggtitle(gene_name) 
    
    #remove_modal_spinner() # remove it when done
  })
  
  
  plotCountRead <- reactive({
    req(filteredintron())
    sum <- filteredintron() %>% 
      group_by(origin) %>%
      summarise(n_read=n())

    return(max(sum$n_read))
  })
  
  plotHeight <- reactive(10*plotCountRead()) 
  
  output$plotreads.ui <- renderUI({
    plotOutput("plot_reads", height = plotHeight())
  })
  
  output$plot_polyA_site <- renderPlot({
    req(AGI_data())
    coords_df <- AGI_data()$COORDS_DF %>% 
      group_by(polya_start_base, polya_length, origin) %>%
      summarise(n=n())
    #ggplot(coords_df, aes(x=polya_start_base, y=round(polya_length), size=n)) +
    #  geom_point(aes(colour = factor(origin)))

    
    ggplot(coords_df, aes(round(polya_length), polya_start_base)) + 
      geom_tile(aes(fill = n)) + 
      geom_text(aes(label = n), 
                size = 3) + 
      coord_fixed() + 
      scale_fill_viridis_c() + 
      guides(fill = FALSE) +
      scale_color_viridis_c()
    
  })
  
  
  ### --- OUTPUTS
  
  

} # /server

#######################################################################
### App calling       
shinyApp(ui, server)

