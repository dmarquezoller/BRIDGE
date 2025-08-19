#' @export
server_function <- function(input, output, session, db_path){

  #### DATA RETRIEVING SERVER ####

  con <- connect_db(db_path)
  session$onSessionEnded(function() {
    try(DBI::dbDisconnect(con), silent = TRUE)
  })
  
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

  
  # Rendering of which table to choose in function of species selected
  output$table_selector <- shiny::renderUI({  
    tables <- valid_tables()
    if (is.null(tables)) {
      return(shiny::helpText("No tables found for this species."))
    } else {
      shiny::selectInput("selected_table", "Select a data table:", choices = tables)
    }
  })
  
  # obtain desired metadata table from the species and table selected 
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
    annotation_table <- input$annotation_table
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

  valid_anno <- shiny::reactive({
    shiny::req(input$species)
    species <- tolower(input$species)
    
    # Query metadata table
    metadata <- DBI::dbGetQuery(con, "SELECT table_name FROM annotation_metadata")
    matches <- grep(paste0("^", species, "_"), metadata$table_name, value = TRUE)
    
    if (length(matches) == 0) {
      return(NULL)
    } else {
      return(matches)
    }
  })


  output$annotation_selector <- shiny::renderUI({
    anno_files <- valid_anno()
    if (is.null(anno_files)) {
      return(shiny::helpText("No annotation files found for this species."))
    } else {
      shiny::selectInput("annotation_table", "Select an annotation file:", choices = anno_files)
    }
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
      rv$species[[table_id]] <- input$species
      rv$annotation[[table_id]] <- annotation_data
      
    }

    cache_key <- paste(table_id, paste(rv$time_cols[[table_id]], collapse = "_"), "dep", sep = "_")
    print(paste("Running analysis for:", table_id, "datatype:", rv$datatype[[table_id]]))

      # Recompute heatmap if the button is clicked

    if (cache$exists(cache_key)) {
        message("Loading DEP output from cache: ", cache_key)
        dep_output <- cache$get(cache_key)
    } else {
        message("Computing and caching DEP output: ", cache_key)
        if (rv$datatype[[table_id]] == 'proteomics') {
          dep_output <- dep2_proteomics(rv$tables[[table_id]], table_id, rv)
        } else if (rv$datatype[[table_id]] == 'phosphoproteomics') {
          dep_output <- dep2_phosphoproteomics(rv$tables[[table_id]], table_id, rv)
        } else if (rv$datatype[[table_id]] == 'rnaseq') {
          dep_output <- dep2_rnaseq(rv$tables[[table_id]], table_id, rv)
        }
        cache$set(cache_key, dep_output)
    }
    print(dep_output)
    rd_names <- colnames( SummarizedExperiment::rowData(dep_output))
    print(rd_names)
    sig_cols <- grep("_significant$", rd_names, value = TRUE)
    valid_contrasts <- sub("_significant$", "", sig_cols) 
    rv$dep_output[[table_id]] <- dep_output
    rv$contrasts[[table_id]] <- valid_contrasts
       
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

  output$all_tables_ui <- renderUI({
  req(length(rv$tables) > 0, !is.null(rv$contrasts))
  tagList(lapply(rv$table_names, function(tbl_name) {
    shinydashboard::box(
      title = paste("Table:", tbl_name), width = 12, solidHeader = TRUE,
      status = "primary", style = "overflow-x: auto", collapsible = TRUE,

      shinydashboard::box(
        title = "Table", width = 12, solidHeader = TRUE, status = "primary",
        style = "overflow-x: auto", collapsible = TRUE,
        h3(), DT::DTOutput(paste0("table_", tbl_name)), h3()
      ),

      shinydashboard::tabBox(
        title = shinyWidgets::actionBttn(paste0("individual_help_", tbl_name),
                 "Help", color = "primary", icon = icon("question-circle"),
                 size = "sm", style = "bordered"),
        width = 12, id = paste0("plotTabs_", tbl_name), selected = "Raw Heatmap",

        tabPanel("Raw Heatmap",
          fluidRow(RawHeatmapUI(paste0("RawHeatmap_", tbl_name), tbl_name))
        ),
        tabPanel("DEP Heatmap",
          fluidRow(DepHeatmapUI(paste0("DepHeatmap_", tbl_name), tbl_name))
        ),
        tabPanel("Volcano Plot",
          fluidRow(VolcanoUI(paste0("Volcano_", tbl_name), tbl_name, rv$contrasts[[tbl_name]]))
        ),
        tabPanel("Gene Expression",
          fluidRow(TimelineUI(paste0("Timeline_", tbl_name), tbl_name))
          # remove TimelineTableUI if merged
        ),
        tabPanel("Enrichment analysis",
          fluidRow(EnrichmentUI(paste0("Enrichment_", tbl_name), tbl_name))
        ),
        tabPanel("PCA",
          fluidRow(pcaUI(paste0("pca_", tbl_name), tbl_name))
        )
      )
    )
  }))
})

  #### END OF MAIN UI FOR EACH TABLE ####
  

  # RENDER OF RAW DATA TABLES
  
  shiny::observe({
    lapply(rv$table_names, function(tbl_name) {
      output[[paste0("table_", tbl_name)]] <- DT::renderDT({
        DT::datatable(rv$tables[[tbl_name]], extensions = "Buttons", options = list(scrollX = TRUE, pageLength = 10, dom = "Bfrtip", buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
      })
    })
  })

 #RawHeatmapServer(paste0("RawHeatmap_", tbl_name), rv) #Function in heatmap_server.R
wired <- reactiveValues(
  timeline = character(),
  raw      = character(),
  dep      = character(),
  pca      = character(),
  enrich   = character()
)

# Helper predicates
has_dep   <- function(tbl) !is.null(rv$dep_output[[tbl]])
has_contr <- function(tbl) !is.null(rv$contrasts[[tbl]]) && length(rv$contrasts[[tbl]]) > 0

observeEvent(list(rv$table_names, rv$dep_output, rv$contrasts), ignoreInit = FALSE, {
  tbls <- rv$table_names

  for (tbl in tbls) {
    # Always-ok modules (need raw table + time cols)
    if (!(tbl %in% wired$timeline) && !is.null(rv$tables[[tbl]]) && !is.null(rv$time_cols[[tbl]])) {
      TimelineServer(paste0("Timeline_", tbl), rv, tbl)
      wired$timeline <- union(wired$timeline, tbl)
    }
    if (!(tbl %in% wired$raw) && !is.null(rv$tables[[tbl]])) {
      RawHeatmapServer(paste0("RawHeatmap_", tbl), rv, tbl)
      wired$raw <- union(wired$raw, tbl)
    }

    # DEP / PCA depend on dep_output
    if (!(tbl %in% wired$dep) && has_dep(tbl)) {
      DepHeatmapServer(paste0("DepHeatmap_", tbl), rv, cache, tbl)
      VolcanoServer(paste0("Volcano_", tbl), rv, tbl)
      wired$dep <- union(wired$dep, tbl)
    }
    if (!(tbl %in% wired$pca) && has_dep(tbl)) {
      pcaServer(paste0("pca_", tbl), rv, tbl)
      wired$pca <- union(wired$pca, tbl)
    }

    # Enrichment also needs contrasts
    if (!(tbl %in% wired$enrich) && has_dep(tbl) && has_contr(tbl)) {
      EnrichmentServer(paste0("Enrichment_", tbl), rv, tbl)
      wired$enrich <- union(wired$enrich, tbl)
    }
  }
})

# Optional cleanup when a table is removed
observeEvent(input$delete_data, {
  wired$timeline <- setdiff(wired$timeline, input$remove)
  wired$raw      <- setdiff(wired$raw, input$remove)
  wired$dep      <- setdiff(wired$dep, input$remove)
  wired$pca      <- setdiff(wired$pca, input$remove)
  wired$enrich   <- setdiff(wired$enrich, input$remove)
})

############################

 integration_ui(input, output, session, rv)

###########################

 help_buttons(input, output, session)

}
