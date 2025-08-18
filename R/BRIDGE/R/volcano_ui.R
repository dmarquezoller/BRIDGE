#' @export
VolcanoUI <- function(id, tbl_name, contrasts) {
  ns <- NS(id)

  # be robust if contrasts is NULL/empty
  if (is.null(contrasts)) contrasts <- character(0)
  sel <- if (length(contrasts)) contrasts[[1]] else NULL

  shinydashboard::box(
    title = "Volcano", width = 12, solidHeader = TRUE, status = "info",
    # SETTINGS
    shiny::fluidRow(
      shinydashboard::box(
        title = "Settings", width = 5, solidHeader = TRUE, status = "info",
        shinyWidgets::virtualSelectInput(
          ns("comparison_volcano"), "Select Comparison:",
          choices = contrasts, selected = sel
        ),
        shiny::selectizeInput(
          ns("volcano_search"), "Search for Gene:", choices = NULL, multiple = TRUE
        ),
        # shinyWidgets::chooseSliderSkin("Flat", color = "#3d8dbc"),
        # −log10(FDR) scale to match pcut <- 10^(-pcut)
        #shiny::sliderInput(
        shiny::numericInput(
          ns("volcano_pcutoff"), "FDR Threshold:",
          min = 0, max = 1, value = 0.05, step = 0.01
        ),
        shiny::numericInput(
          ns("volcano_fccutoff"), "FC Threshold:",
          min = 0, max = 10, value = 1, step = 0.1
        ),
        shinyWidgets::actionBttn(
          ns("compute_volcano"),
          shiny::span("Compute Volcano", style = "color: white;"),
          style = "simple", color = "primary", size = "sm"
        )
      ),
      shinydashboard::box(
        title = "Volcano Plot", width = 7, solidHeader = TRUE, status = "info",
        uiOutput(ns("volcano_slot"))
      )
    ),
    h5(),
    shinydashboard::box(
      title = "Significant Genes", width = 12, solidHeader = TRUE,
      status = "info", collapsible = TRUE, collapsed = FALSE,
      DT::DTOutput(ns("volcano_sig_table"))
    )
  )
}