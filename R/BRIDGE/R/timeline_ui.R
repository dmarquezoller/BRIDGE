
#' @export
timeline_gene_search <- function(tbl_name, rv) {
    if (rv$datatype[[tbl_name]] == "rnaseq") {
        scale_choices <- c("Continous", "Log-scale", "TPM", "FPKM", "TMM", "CPM")
    } else {
        scale_choices <- c("Continous", "Log-scale", "Total Intensity", "Median Normalization")
    }
        shinydashboard::box(
        title = "Gene Search", width = 6, solidHeader = TRUE, status = "info", 
        selectizeInput(
            inputId = paste0("search_gene_", tbl_name),
            label = "Search your gene of interest:",
            choices = NULL,
            multiple = TRUE
        ),
        shiny::selectInput(
            inputId = paste0("scale_", tbl_name),
            label = "Select scale:",
            choices = scale_choices,
            selected = "Continous"
        )
    )
}

#' @export
timeline_plot <- function(tbl_name) {
    shinydashboard::box(
        title="Expression Timeline", width = 12, solidHeader = T, status = "info",   
        h3(),
        plotOutput(paste0("time_plot_", tbl_name)),
        h3()
    )
}

#' @export
timeline_table <- function(tbl_name) {
    shinydashboard::box(
        title="Table Entries", width = 12, solidHeader = T, status = "info", style= "overflow-x: auto", collapsible = T, collapsed = F, 
        h3(), 
        DTOutput(paste0("time_plot_dt_", tbl_name), height = "300px"), 
        h3()
    )
}