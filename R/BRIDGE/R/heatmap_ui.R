#' @export
heatmap_ui <- function(id) {
  ns <- NS(id)
  shinydashboard::tabBox(
    title = "Heatmaps",
    shiny::tabPanel("Raw Heatmap", plotOutput(ns("raw_ht"))),
    shiny::tabPanel("DEP Heatmap", plotOutput(ns("dep_ht")))
  )
}

#' @export
raw_heatmap_ui <-function(tbl_name) { 
                    shinydashboard::box(
                      title = "Raw Data Heatmap", width = 12, solidHeader = T, status = "info", 
                        shiny::fluidRow(
                            shiny::column(
                              width = 7, 
                              shinycssloaders::withSpinner(
                                plotOutput(paste0("raw_ht_", tbl_name), width = "500px", height = "500px"), 
                                type = 8, 
                                color = "#2b8cbe",
                                caption = "Loading..."
                              )
                            ),
                            shiny::column(
                              width = 5, 
                              h1(), 
                              shinyWidgets::actionBttn(paste0("compute_raw_ht_", tbl_name), shiny::span("Compute Raw Heatmap", style = "color: white;"), style = "simple", color = "primary", size = "sm"),
                              h3(),
                              p("Have in mind that for big heatmaps it may take long.")
                            )
                          )
                        )
                  }

#' @export
dep_heatmap_ui  <- function(tbl_name) { 
                        shinydashboard::box(
                         title = "Heatmap", width = 12, solidHeader = TRUE, status = "info",
                         shiny::fluidRow(),
                         shiny::fluidRow(
                           shiny::column(
                             width = 7,
                             shinycssloaders::withSpinner(
                              plotOutput(paste0("ht_", tbl_name), width = "500px", height = "500px"),
                              type = 8,
                              color = "#2b8cbe", 
                              caption = "Loading..."
                             )
                           ),
                           shiny::column(
                             width = 5,
                             h1(),
                             p("Enable clustering:"),
                             shinyWidgets::switchInput(paste0("clustering_", tbl_name), value = T, onLabel = "YES", offLabel = "NO", width = 'auto'),
                             shiny::numericInput(paste0("num_clusters_", tbl_name), "Select Number of Clusters", min = 2, step = 1, value = 3),
                             shinyWidgets::actionBttn(paste0("recompute_heatmap_", tbl_name), shiny::span("Compute Heatmap", style = "color: white;"), style = "simple", color = "primary", size = "sm"),
                             h1(),
                             uiOutput(paste0("optimal_k", tbl_name))                            
                           )
                         ),
                         h5(),
                         shinydashboard::box(title="Selected Genes", width = 12, solidHeader = T, status = "info", style= "overflow-x: auto",collapsible = T, collapsed = F,
                             h3(),
                             DT::DTOutput(paste0("ht_sig", tbl_name), height = "300px"),
                             h3()
                             
                         )
                       )
                    }