#' @export
pcaUI <- function(id, tbl_name) {
  ns <- NS(id)
  shinydashboard::box(
    title = "PCA", width = 12, solidHeader = TRUE, status = "primary",
    actionButton(ns("compute"), "Compute PCA"),
    plotOutput(ns("plot"), height = "480px")
  )
}
