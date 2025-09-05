#' @export
data_server <- function(id, con, rv) {
    shiny::moduleServer(id, function(input, output, session) {
        output$all_tables_ui <- shiny::renderUI({
            shiny::req(length(rv$tables) > 0)
            lapply(rv$table_names, function(tbl_name) {
                shinydashboard::box(title = paste("Table:", tbl_name), DT::DTOutput(paste0("table_", tbl_name)))
            })
        })
    })
}
