#' @export
pca_ui <- function(tbl_name) {
  box(
    title = "PCA Analysis", width = 12, solidHeader = TRUE, status = "info",
    fluidRow(
      column(
        width = 12,
        div(
          style = "text-align: center; margin-bottom: 20px;",
          actionBttn(
            inputId = paste0("compute_pca_", tbl_name),
            label = span("Compute PCA", style = "color: white;"),
            style = "simple",
            color = "primary",
            size = "md"
          )
        )
      )
    ),
    fluidRow(
      column(
        width = 12,
        div(
          style = "padding: 10px;",
          plotOutput(outputId = paste0("pca_", tbl_name)),
        )
      )
    )
  )
}
