#' @export
heatmap_ui <- function(id) {
  ns <- NS(id)
  tabBox(
    title = "Heatmaps",
    tabPanel("Raw Heatmap", plotOutput(ns("raw_ht"))),
    tabPanel("DEP Heatmap", plotOutput(ns("dep_ht")))
  )
}

#' @export
raw_heatmap_ui <-function(tbl_name) { 
                    box(
                      title = "Raw Data Heatmap", width = 12, solidHeader = T, status = "info", 
                        fluidRow(
                            column(
                              width = 7, 
                              withSpinner(
                                plotOutput(paste0("raw_ht_", tbl_name), width = "500px", height = "500px"), 
                                type = 8, 
                                color = "#2b8cbe",
                                caption = "Loading..."
                              )
                            ),
                            column(
                              width = 5, 
                              h1(), 
                              actionBttn(paste0("compute_raw_ht_", tbl_name), span("Compute Raw Heatmap", style = "color: white;"), style = "simple", color = "primary", size = "sm"),
                              h3(),
                              p("Have in mind that for big heatmaps it may take long.")
                            )
                          )
                        )
                  }

#' @export
dep_heatmap_ui  <- function(tbl_name) { 
                        box(
                         title = "Heatmap", width = 12, solidHeader = TRUE, status = "info",
                         fluidRow(),
                         fluidRow(
                           column(
                             width = 7,
                             withSpinner(
                              plotOutput(paste0("ht_", tbl_name), width = "500px", height = "500px"),
                              type = 8,
                              color = "#2b8cbe", 
                              caption = "Loading..."
                             )
                           ),
                           column(
                             width = 5,
                             h1(),
                             p("Enable clustering:"),
                             switchInput(paste0("clustering_", tbl_name), value = T, onLabel = "YES", offLabel = "NO", width = 'auto'),
                             numericInputIcon(paste0("num_clusters_", tbl_name), "Select Number of Clusters", min = 2, step = 1, value = 3),
                             actionBttn(paste0("recompute_heatmap_", tbl_name), span("Compute Heatmap", style = "color: white;"), style = "simple", color = "primary", size = "sm")
                           )
                         ),
                         h5(),
                         box(title="Selected Genes", width = 12, solidHeader = T, status = "info", style= "overflow-x: auto",collapsible = T, collapsed = F,
                             h3(),
                             DTOutput(paste0("ht_sig", tbl_name), height = "300px"),
                             h3()
                             
                         )
                       )
                    }