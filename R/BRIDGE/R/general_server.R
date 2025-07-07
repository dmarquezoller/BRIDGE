#' @export
server_function <- function(input, output, session, db_path){

  #### DATA RETRIEVING SERVER ####

  con <- connect_db(db_path)
  all_tables <- DBI::dbListTables(con)
  
  #creation of cache starr object
  cache <- storr::storr_dbi(
    tbl_data = "storr_data",
    tbl_keys = "storr_keys",
    con = con
  )
  rv <- reactiveValues(tables=list(), table_names=character(), ht_matrix = list(), time_cols = list(), datatype = character(), constrasts = list()) #variable that stores most of the important values for each table
  
  valid_tables <- shiny::reactive({
    shiny::req(input$species)
    species <- tolower(input$species)
    
    # Query metadata table
    metadata <- DBI::dbGetQuery(con, "SELECT table_name FROM table_metadata")
    matches <- grep(paste0("^", species, "_"), metadata$table_name, value = TRUE)
    
    if (length(matches) == 0) {
      return(NULL)
    } else {
      return(matches)
    }
  })

  
  # Rendering of which table to choose in function of specie selected
  output$table_selector <- shiny::renderUI({  
    tables <- valid_tables()
    if (is.null(tables)) {
      return(shiny::helpText("No tables found for this species."))
    } else {
      shiny::selectInput("selected_table", "Select a data table:", choices = tables)
    }
  })
  
  # obtain desired metadata table from the specie and table selected 
  table_metadata <- shiny::reactive({
    shiny::req(input$selected_table)
    query <- sprintf("SELECT * FROM table_metadata WHERE table_name = '%s'", input$selected_table)
    metadata <- DBI::dbGetQuery(con, query)
    
    if (nrow(metadata) == 0) return(NULL)
    
    # Get all non-null fields in the metadata row
    column_types <- names(metadata)[!is.na(metadata[1, ]) & metadata[1, ] != ""]
    
    # Exclude 'table_name'
    column_types <- setdiff(column_types, "table_name")
    
    # Create a named list with each type and its corresponding columns
    meta_list <- lapply(column_types, function(ct) unlist(strsplit(metadata[[ct]], ",")))
    names(meta_list) <- column_types
    
    return(meta_list)
  })

  annotation_metadata <- shiny::reactive({
    shiny::req(input$species)
    annotation_table <- paste0(input$species, "_annotation")
    query <- sprintf("SELECT * FROM '%s'", annotation_table)
    return(query)
  })
  
  # COLUMN SELECTOR INPUT 

  # From the metadata table retrieved show the possible timepoints columns to choose
  output$column_selector <-shiny::renderUI({
    meta <- table_metadata()
    if (is.null(meta)) return(shiny::helpText("No metadata available for this table."))
    
    timepoint_cols <- unlist(strsplit(meta$timepoint_columns, ","))
    
    shinyWidgets::pickerInput(
      inputId = "timepoints_selected",
      label = "Select timepoint columns to load:",
      choices = c(timepoint_cols),
      selected = NULL,
      multiple = TRUE,
      options = shinyWidgets::pickerOptions(
        actionsBox = TRUE,
        liveSearch = TRUE,
        noneSelectedText = "Select columns",
        deselectAllText = "None",
        selectAllText = "Select All",
        dropupAuto = FALSE
      ),
      width = "100%"
    )
  })
  
  # LOAD DATA BUTTON SERVER
  
  # Server side of loading the data, selecting the desired columns 
  shiny::observeEvent(input$load_data, {
    shiny::req(input$selected_table, input$timepoints_selected)

    meta <- table_metadata()
    annotation <- annotation_metadata()
    if (is.null(meta)) return()
    # Always include identifiers
    id_cols <- trimws(unlist(strsplit(meta$identifier_columns, ",")))
    timepoint_cols <- trimws(input$timepoints_selected)
    
    safe_cols <- sprintf('"%s"', unique(c(id_cols, timepoint_cols)))
    col_string <- paste(safe_cols, collapse = ", ")
    query <- sprintf("SELECT %s FROM %s", col_string, input$selected_table)
    new_data <- DBI::dbGetQuery(con, query, )
    annotation_data <- DBI::dbGetQuery(con, annotation, )
    table_id <- input$selected_table
    
    datatype <- strsplit(input$selected_table, "_")[[1]][2]
    if (!(table_id %in% rv$table_names)) {
      rv$tables[[table_id]] <- new_data
      rv$table_names <- c(rv$table_names, table_id)
      matrix_data <- new_data[, timepoint_cols]
      rv$ht_matrix[[table_id]] <- as.matrix(matrix_data)
      rv$id_cols[[table_id]] <- id_cols
      rv$time_cols[[table_id]] <- input$timepoints_selected
      rv$datatype[[table_id]] <- datatype
      rv$annotation[[table_id]] <- annotation_data
      rv$species[[table_id]] <- input$species
      
    }
  })
  
  # DELETE DATA SELECTOR

  shiny::observe({
    shiny::updateSelectInput(session, "remove", choices = rv$table_names, selected = NULL)
  })
  
  # DELETE DATA BUTTON SERVER

  shiny::observeEvent(input$delete_data, {
    shiny::req(input$remove)
    
    if (input$remove %in% rv$table_names) {
      # Remove the table from reactive values
      rv$tables[[input$remove]] <- NULL
      rv$table_names <- setdiff(rv$table_names, input$remove)
      
      # Optionally reset the dropdown selection
      shiny::updateSelectInput(session, "remove", choices = rv$table_names, selected = NULL)
    }
  })
  
  # LOADED TABLES INFO

  output$loaded_info <- shiny::renderUI({
    if (length(rv$table_names) == 0) {
      return("No table loaded yet.")
    } else {
      # Create an unordered list with each table name as a list item
       htmltools::tagList(
        shiny::tags$ul(
          lapply(rv$table_names, function(tbl_name) {
            shiny::tags$li(tbl_name)
          })
        )
      )
    }
  })
  
  #### END OF DATA RETRIEVING SERVER ####
  

  #### MAIN UI FOR EACH TABLE #### (includes functions for heatmap ui)

  output$all_tables_ui <- shiny::renderUI({
    shiny::req(length(rv$tables) > 0)

    lapply(rv$table_names, function(tbl_name) {
      shinydashboard::box(title=paste("Table: ", tbl_name), width=12, solidHeader = TRUE, status="primary", style="overflow-x: auto", collapsible = T, collapsed = F,
          shinydashboard::box(title="Table", width = 12, solidHeader = T, status = "primary", style= "overflow-x: auto", collapsible = T, collapsed = F,
              h3(),
              DT::DTOutput(paste0("table_", tbl_name)),
              h3()              
          ),
          shinydashboard::tabBox(
            title = shinyWidgets::actionBttn("individual_help", "Help", color = "primary" ,icon=shiny::icon("question-circle"), size="sm", style = "bordered"), width = 12,
            id = "plotTabs", selected = "Raw Heatmap",
            shiny::tabPanel("Raw Heatmap", 
                     shiny::fluidRow(
                        raw_heatmap_ui(tbl_name) #Function in heatmap_ui.R
                      )
                    
            ),
            shiny::tabPanel("DEP Heatmap",
                     shiny::fluidRow(
                        dep_heatmap_ui(tbl_name) #Function in heatmap_ui.R
                     )
            ),
            shiny::tabPanel("Volcano Plot",
                     shiny::fluidRow(
                       shinydashboard::box(title = "Settings", width = 5, solidHeader = T, status = "info",
                           shinyWidgets::virtualSelectInput(paste0("comparison_volcano_", tbl_name), "Select Comparison:", 
                                              choices = rv$contrasts[[tbl_name]], 
                                              selected = rv$contrasts[[tbl_name]][1]
                                              ),
                           shiny::selectizeInput(paste0("volcano_search_",tbl_name), "Search for Gene:",  choices = NULL, multiple=T),
                           shinyWidgets::chooseSliderSkin("Flat", color="#3d8dbc"),
                           shiny::numericInput(paste0("volcano_pcutoff_", tbl_name), "FDR Threshold:", 
                                       min = 0, max = 10, value = 1.3, step = 0.1),
                           shiny::numericInput(paste0("volcano_fccutoff_", tbl_name), "FC Threshold:", 
                                       min = 0, max = 10, value = 1, step = 0.1),
                           shinyWidgets::actionBttn(paste0("compute_volcano_", tbl_name), shiny::span("Compute  Volcano", style = "color: white;"), style = "simple", color = "primary", size = "sm")
                       ),
                       shinydashboard::box(
                         title = "Volcano Plot", width = 7, solidHeader = TRUE, status = "info",
                         shinycssloaders::withSpinner(
                          plotly::plotlyOutput(paste0("volcano_", tbl_name)),
                          type = 8,
                          color = "#2b8cbe", 
                          caption = "Loading..."
                         )
                       )
                       
                     ),
                     shinydashboard::box(title="Significant Genes", width=12, solidHeader = T, status="info", collapsible = T, collapsed = F,
                         h3(),
                         DT::DTOutput(paste0("volcano_sig_table_", tbl_name)),
                         h3()
                     )
            ),
            shiny::tabPanel("Gene Expression", 
                     shiny::fluidRow(
                        timeline_gene_search(tbl_name, rv) #Function for the ui of the timeline search in timeline_ui.R
                     ),
                     shiny::fluidRow(
                        timeline_plot(tbl_name) #Function for the ui of the timeline plot in timeline_ui.R
                     ), 
                     shiny::fluidRow(
                        timeline_table(tbl_name) #Function for the ui of the timeline table in timeline_ui.R
                     )
            ),
            shiny::tabPanel("Enrichment analysis", 
                     shiny::fluidRow(
                        enrichment_settings_ui(tbl_name, rv), 
                        enrichment_plots_ui(tbl_name)
                     )
            ),
            shiny::tabPanel("PCA", 
                    shiny::fluidRow(
                      pca_ui(tbl_name)
                    )
            )
          )
      )
      
    })
  })

  #### END OF MAIN UI FOR EACH TABLE ####
  

  # RENDER OF RAW DATA TABLES
  
  shiny::observe({
    lapply(rv$table_names, function(tbl_name) {
      output[[paste0("table_", tbl_name)]] <- DT::renderDT({
        DT::datatable(rv$tables[[tbl_name]],options = list(scrollX = TRUE, pageLength = 10))
      })
    })
  })



####### RAW HEATMAP ########

 raw_heatmap_server(input, output, session, rv) #Function in heatmap_server.R

####### DEP HEATMAP ########
 
 dep_heatmap_server(input, output, session, rv, cache) #Function in heatmap_server.R

###### TIMELINE PLOT #######

 timeline_server(input, output, session, rv) #Function in timeline_server.R
  
######### PCA PLOT #########

 pca_server(input, output, session, rv) #Function in pca_server.R

##### ENRICHMENT PLOT ######

 enrichment_server(input, output, session, rv)

############################

 integration_ui(input, output, session, rv)

###########################

 help_buttons(input, output, session)

}
