#' @export
connect_db <- function(db_path) {
  dbConnect(SQLite(), db_path)
}

#' @export
list_tables <- function(con) {
  dbListTables(con)
}