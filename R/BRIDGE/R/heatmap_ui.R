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
RawHeatmapUI <- function(id, tbl_name) {
  ns <- NS(id)
  shinydashboard::box(
    title = "Raw Data Heatmap", width = 12, solidHeader = TRUE, status = "info",
    fluidRow(
      column(
        5,
        shinyWidgets::actionBttn(
          ns("compute"),
          span("Compute Raw Heatmap", style = "color:white;"),
          style = "simple", color = "primary", size = "sm"
        ),
        p("Might take long for first time computation.")
      ),
      column(
        7,
        # IMPORTANT: don't wrap uiOutput in a spinner here
        uiOutput(ns("plot_slot"))
      )
    )
  )
}


#' @export
DepHeatmapUI <- function(id, tbl_name) {
  ns <- NS(id)
  shinydashboard::box(
    title = "Heatmap", width = 12, solidHeader = TRUE, status = "info",
    fluidRow(
      column(
        width = 7,
        uiOutput(ns("ht_slot"))
      ),
      column(
        width = 5,
        h1(), p("Change thresholds"),
        # Use −log10(FDR) to match p_cut <- 10^(−input$heatmap_pcutoff)
        #sliderInput(ns("heatmap_pcutoff"), "FDR Threshold (−log10):",
        numericInput(ns("heatmap_pcutoff"), "FDR Threshold:",
                    min = 0, max = 1, value = 0.05, step = 0.01),
        numericInput(ns("heatmap_fccutoff"), "FC Threshold:",
                     min = 0, max = 10, value = 1, step = 0.1),
        p("Enable clustering:"),
        shinyWidgets::switchInput(ns("clustering"), value = TRUE,
                                  onLabel = "YES", offLabel = "NO", width = 'auto'),
        numericInput(ns("num_clusters"), "Select Number of Clusters",
                     min = 2, step = 1, value = 3),
        shinyWidgets::actionBttn(
          ns("recompute_heatmap"),
          span("Compute Heatmap", style = "color: white;"),
          style = "simple", color = "primary", size = "sm"
        ),
        h1(),
        uiOutput(ns("optimal_k"))   # <- fixed id (no paste0 with tbl_name)
      )
    ),
    h5(),
    shinydashboard::box(
      title = "Selected Genes", width = 12, solidHeader = TRUE,
      status = "info", style = "overflow-x: auto",
      collapsible = TRUE, collapsed = FALSE,
      DT::DTOutput(ns("ht_sig"), height = "300px")   # <- fixed id
    )
  )
}