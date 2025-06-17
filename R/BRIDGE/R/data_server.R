#' @export
data_server <- function(id, con, rv) {
  moduleServer(id, function(input, output, session) {
    output$all_tables_ui <- renderUI({
      req(length(rv$tables) > 0)
      lapply(rv$table_names, function(tbl_name) {
        box(title = paste("Table:", tbl_name), DTOutput(paste0("table_", tbl_name)))
      })
    })
  })
}