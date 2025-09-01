#' @export
pcaUI <- function(id, tbl_name) {
    ns <- NS(id)

    shinydashboard::box(
        title = "PCA", width = 12, solidHeader = TRUE, status = "info",
        shinyWidgets::actionBttn(ns("compute"),
            span("Compute PCA", style = "color: white;"),
            style = "simple", color = "primary", size = "sm"
        ),
        shinycssloaders::withSpinner(
            plotOutput(ns("plot"), height = "480px"),
            type = 8, color = "#2b8cbe", caption = "Loading..."
        ),
        h5(),
        shinydashboard::box(
            title = "PCA loadings", width = 12, solidHeader = TRUE,
            status = "info", style = "overflow-x: auto",
            collapsible = TRUE, collapsed = FALSE,
            DT::DTOutput(ns("pcs"), height = "300px") 
        )
    )
}
