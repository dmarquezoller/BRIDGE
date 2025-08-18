#' @export
TimelineGeneSearchUI <- function(id, tbl_name, rv) {
    ns <- NS(id)
    if (rv$datatype[[tbl_name]] == "rnaseq") {
        scale_choices <- c("Continous", "Log-scale", "TPM", "FPKM", "TMM", "CPM")
    } else {
        scale_choices <- c("Continous", "Log-scale", "Total Intensity", "Median Normalization")
    }
        shinydashboard::box(
        title = "Gene Search", width = 6, solidHeader = TRUE, status = "info", 
        selectizeInput(
            inputId = ns(paste0("search_gene_", tbl_name)),
            label = "Search your gene of interest:",
            choices = NULL,
            multiple = TRUE
        ),
        shiny::selectInput(
            inputId = ns(paste0("scale_", tbl_name)),
            label = "Select scale:",
            choices = scale_choices,
            selected = "Continous"
        )
    )
}

#' @export
TimelinePlotUI <- function(id, tbl_name) {
    ns <- NS(id)
    shinydashboard::box(
        title="Expression Timeline", width = 12, solidHeader = T, status = "info",   
        h3(),
        plotOutput(ns(paste0("time_plot_", tbl_name))),
        h3()
    )
}

#' @export
TimelineTableUI <- function(id, tbl_name) {
    ns <- NS(id)
    shinydashboard::box(
        title="Table Entries", width = 12, solidHeader = T, status = "info", style= "overflow-x: auto", collapsible = T, collapsed = F, 
        h3(), 
        DT::DTOutput(ns(paste0("time_plot_dt_", tbl_name)), height = "300px"), 
        h3()
    )
}