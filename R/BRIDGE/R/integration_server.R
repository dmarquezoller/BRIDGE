#source("R/global.R")
#source("R/modules/server/int_timeline_server.R")
#source("R/modules/server/int_heatmap_server.R")
#source("R/modules/helpers/integration_helper.R")
#source("R/modules/server/raw_integration_server.R")
#source("R/modules/server/processed_integration_server.R")

#' @export
integration_ui <- function(input, output, session, rv) {
  combined_data <- reactiveVal(NULL)

  output$integration_ui <- renderUI({
    req(length(rv$tables) > 0)
    tagList(
        tabBox(title = "Integration", id = "integration_tabs", selected = "Raw Integration", width = 12, 
            tabPanel("Raw Integration",
                fluidRow(
                    box(title = "Integration Settings", width = 12, status = "success", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                    fluidRow(
                        box(
                        title = "Integrate", width = 4, solidHeader = TRUE, status = "primary",
                        tags$div(style = "height:8px;"),
                        div(style = "text-align: center;",    
                            actionBttn("raw_int_help", "Help", color = "primary" ,icon=icon("question-circle"), size="sm", style = "bordered")
                        ),
                        tags$div(style = "height:8px;"),
                        pickerInput(inputId = "integration", label = "Integrate", choices = rv$table_names, multiple = TRUE, options = pickerOptions(container = "body"), width = "100%"),
                        uiOutput("integration_col_selector"),
                        div(style = "text-align: center;",
                            actionBttn("integrate_data", span("Integrate", style = "color: black;"), icon = span(icon("arrow-right-to-bracket"), style = "color: black;"), style = "jelly", size = "sm", class = "btn-primary")
                        )
                        ),
                        uiOutput("preview_box")
                    )
                ), 
                box(
                    title = "Integrated Raw Table", width = 12, solidHeader = TRUE, status = "info", collapsible = TRUE, collapsed = FALSE, 
                    DTOutput("integration_combined_table")
                ),
                box(
                    title = "Integrated Timeline Plot", width = 12, solidHeader = TRUE, status = "info",
                    selectizeInput(
                        inputId = "search_gene_integration",
                        label = "Search your gene of interest:",
                        choices = NULL,
                        multiple = TRUE
                    ),
                    selectInput(
                        inputId = "scale_integration",
                        label = "Select scale:",
                        choices = c("Continous", "Log-scale"),
                        selected = "Continous"
                    ),
                    plotOutput("integration_timeline_plot")
                )
            )
            ),
            tabPanel("Processed Integration",
                fluidRow(
                    box(title = "Processed Integration Settings", width = 12, status = "success", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                    fluidRow(
                        box(
                            title = "Process Integrate", width = 4, solidHeader = TRUE, status = "primary",
                            tags$div(style = "height:8px;"),
                            div(style = "text-align: center;",
                                actionBttn("processed_int_help", "Help", color = "primary" ,icon=icon("question-circle"), size="sm", style = "bordered")
                            ),
                            tags$div(style = "height:8px;"),
                            pickerInput(
                                inputId = "processed_integration", label = "Integrate", 
                                choices = rv$table_names, multiple = TRUE, 
                                options = pickerOptions(container = "body"), width = "100%"
                            ),
                            p("Select for each table upon which contrast do the filtering (they should match)"),
                            uiOutput("comparison_col_selector_pi"),
                            p("----------"),
                            numericInput("heatmap_k", "Number of clusters (k):", value = 3, min = 2, max = 10),
                            numericInput("pval_thresh_pi", "P-value threshold:", value = 0.05, min = 0, step = 0.01),
                            numericInput("lfc_thresh_pi", "LFC threshold:", value = 1, min = 0, step = 0.1),

                            div(style = "text-align: center;",
                                actionBttn("process_integrate_data", span("Integrate", style = "color: black;"), 
                                            icon = span(icon("arrow-right-to-bracket"), style = "color: black;"), 
                                            style = "jelly", size = "sm", class = "btn-primary")
                            )
                            ),
                            uiOutput("preview_box_integrate")
                    )),
                    box(title = "Tables", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, 
                        fluidRow(
                            uiOutput("processed_integrated_tables")
                        )
                    ),
                    box(title = "Integrated Heatmaps", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                      uiOutput("integrated_heatmaps")
                    ),
                    box(title = "LFC Scatterplot", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                      uiOutput("lfc_scatter_ui")
                    )

                ))
       )
    )
  })

  # Preview box with column matching and unified table
  output$preview_box <- renderUI({
    req(length(input$integration) > 0)
    tagList(
      box(
        title = "Integration Preview", width = 8, solidHeader = TRUE, status = "primary",
        uiOutput("integration_column_matching")
      )
    )
  })



  output$comparison_col_selector_pi <- renderUI({
    req(input$processed_integration)
    lapply(input$processed_integration, function(tbl) {
      available_comparisons <- rv$contrasts[[tbl]]
      selectInput(
        inputId = paste0("pi_comparison_selected_", tbl),
        label = paste0("Select comparison for ", tbl),
        choices = available_comparisons,
        selected = NULL,
        width = "100%", 
        multiple = FALSE
      )
    })
  })



    output$preview_box_integrate <- renderUI({
      req(rv$integration_preview_dims)
      req(rv$optimal_k)

      dims_ui <- lapply(names(rv$integration_preview_dims), function(tbl) {
        dims <- rv$integration_preview_dims[[tbl]]
        tagList(
          tags$li(tags$b(tbl)),
          tags$ul(
            tags$li(paste("Original:", paste(dims$original, collapse = " x "))),
            tags$li(paste("After filtering:", paste(dims$filtered, collapse = " x "))),
            tags$li(paste("After intersection:", paste(dims$intersected, collapse = " x ")))
          )
        )
      })

      tagList(
        box(
          title = "Processed Integration Preview", width = 8, solidHeader = TRUE, status = "primary",
          tags$ul(
            tags$li(tags$b("Table dimensions at each step:")),
            dims_ui,
            tags$li(tags$b("Selected comparisons:")),
            tags$ul(
              lapply(input$processed_integration, function(tbl) {
                comp <- input[[paste0("pi_comparison_selected_", tbl)]]
                tags$li(paste(tbl, ":", comp))
              })
            ),
            tags$li(tags$b("Optimal k following the elbow rule:")),
            tags$ul(
              tags$li(paste("The optimal k is: ", rv$optimal_k))
            )
          )
        )
      )
    })


    output$processed_integrated_tables <- renderUI({
        req(rv$intersected_tables_processed)

        tables <- rv$intersected_tables_processed

        # Generate a UI list of DTOutput placeholders
        tagList(
            lapply(names(tables), function(tbl) {
            box(
                title = paste("Processed Table:", tbl),
                width = 12,
                solidHeader = TRUE,
                status = "info",
                DTOutput(outputId = paste0("processed_tbl_", tbl))
            )
            })
        )
    })

    observe({
        req(rv$intersected_tables_processed)
        lapply(names(rv$intersected_tables_processed), function(tbl) {
            local({
            table_name <- tbl
            output_id <- paste0("processed_tbl_", table_name)
            data <- as.data.frame(rv$intersected_tables_processed[[table_name]])
            output[[output_id]] <- renderDT({
                datatable(data, options = list(scrollX = TRUE, pageLength = 5))
            })
            })
        })
    })




  # Column selector for each selected table
  output$integration_col_selector <- renderUI({
    req(length(input$integration) > 0)
    lapply(input$integration, function(int_tbl) {
      cols <- rv$time_cols[[int_tbl]]
      pickerInput(
        inputId = paste0("cols_selected_int", int_tbl),
        label = paste0("Select columns from ", int_tbl, " to load:"),
        choices = c(cols),
        selected = NULL,
        multiple = TRUE,
        options = pickerOptions(
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
  })




  # Show how columns will be matched
  output$integration_column_matching <- renderUI({
    req(length(input$integration) > 1)

    selected_columns <- lapply(input$integration, function(tbl) {
      input[[paste0("cols_selected_int", tbl)]]
    })

    if (length(unique(sapply(selected_columns, length))) > 1) {
      return(p("Column numbers differ. Matching not possible.", style = "color: red;"))
    }

    n <- length(selected_columns[[1]])
    if (n == 0) return(NULL)

    matches <- lapply(seq_len(n), function(i) {
      paste(
        sapply(seq_along(selected_columns), function(j) {
          selected_columns[[j]][i]
        }),
        collapse = " ⇄ "
      )
    })

    tagList(
      tags$ul(
        lapply(matches, function(line){
          tags$li(line)
        })
      )
    )
  })

    ### RAW INTEGRATION PIPELINE

    raw_integration(input, output, session, rv, combined_data)

    ### PROCESSED INTEGRATION PIPELINE

    processed_integration(input, output, session, rv)

    ### HEATMAP PLOTING PIPELINE

    int_heatmap_server(input, output, session, rv)



  # Render integrated table
    output$integration_combined_table <- renderDT({
    if (is.null(combined_data())) {
        return(
        datatable(
            data.frame(Message = "Click Integrate button to see the Integrated Table"),
            options = list(
            dom = 't',
            ordering = FALSE,
            paging = FALSE,
            searching = FALSE
            ),
            rownames = FALSE
        )
        )
    }

    DT::datatable(combined_data(), options = list(scrollX = TRUE, pageLength = 10))
    })



}
