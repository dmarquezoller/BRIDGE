#' @export
EnrichmentUI <- function(id, tbl_name) {
  ns <- NS(id)
  tagList(
    shinydashboard::box(
      title = "Enrichment settings", width = 12, solidHeader = TRUE, status = "primary",
      fluidRow(
        column(4, selectInput(ns("comparison_db"), "Database", choices = c("GO", "KEGG", "Reactome"))),
        column(4, uiOutput(ns("contrast_ui"))),
        column(4, actionButton(ns("compute_enrichment"), "Compute enrichment"))
      )
    ),
    shinydashboard::box(
      title = "Enrichment results", width = 12, solidHeader = TRUE, status = "primary",
      uiOutput(ns("enrichment"))   # will insert plotOutput(ns("enrichment_plot")) or a message
    )
  )
}