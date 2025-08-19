#' @export
EnrichmentUI <- function(id, tbl_name) {
  ns <- NS(id)
  tagList(
    shinydashboard::box(
      title = "Enrichment settings", width = 12, solidHeader = TRUE, status = "info",
      fluidRow(
        column(4, selectInput(ns("comparison_db"), "Database", choices = c("GO", "KEGG", "Reactome"))),
        column(4, uiOutput(ns("contrast_ui"))),
        column(4, shinyWidgets::actionBttn(ns("compute_enrichment"), span("Compute Enrichment", style = "color: white;"), style = "simple", color = "primary", size = "sm" ))
      )
    ),
    shinydashboard::box(
      title = "Enrichment results", width = 12, solidHeader = TRUE, status = "info",
      uiOutput(ns("enrichment"))   # will insert plotOutput(ns("enrichment_plot")) or a message
    )
  )
}