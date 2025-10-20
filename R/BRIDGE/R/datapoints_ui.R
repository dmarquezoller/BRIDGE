datapointsUI <- function(id, tbl_name) {
    ns <- NS(id)
    tagList(
        shinydashboard::box(
            title = "datapoints", width = 12, solidHeader = TRUE, status = "info",
            fluidRow(
                column(4, selectizeInput(ns("search_gene"), "Gene(s)", choices = NULL, multiple = TRUE)),
                column(4, selectInput(ns("scale"), "Scale",
                    choices = c(
                        "Continous", "Log-scale", "Median Normalization",
                        "Total Intensity", "FPKM", "TPM", "TMM", "CPM"
                    ),
                    selected = "Continous"
                ))
            ),
            hr(),
            # DT::DTOutput(ns("data_plot_dt")),
            # br(),
            plotOutput(ns("data_plot"), height = "520px")
        )
    )
}

#' @export
datapointsTableUI <- function(id, tbl_name) {
    ns <- NS(id)
    shinydashboard::box(
        title = "Table Entries", width = 12, solidHeader = T, status = "info", style = "overflow-x: auto", collapsible = T, collapsed = F,
        h3(),
        # DT::DTOutput(ns(paste0("data_plot_dt_", tbl_name)), height = "300px"),
        DT::DTOutput(ns("data_plot_dt")),
        h3()
    )
}
