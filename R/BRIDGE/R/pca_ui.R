#' @export
pca_ui <- function(tbl_name) {
  shinydashboard::box(
    title = "PCA Analysis", width = 12, solidHeader = TRUE, status = "info",
    shiny::fluidRow(
      shiny::column(
        width = 12,
       shiny::div(
          style = "text-align: center; margin-bottom: 20px;",
          shinyWidgets::actionBttn(
            inputId = paste0("compute_pca_", tbl_name),
            label = shiny::span("Compute PCA", style = "color: white;"),
            style = "simple",
            color = "primary",
            size = "md"
          )
        )
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 12,
       shiny::div(
          style = "padding: 10px;",
          plotOutput(outputId = paste0("pca_", tbl_name)),
        )
      )
    )
  )
}
