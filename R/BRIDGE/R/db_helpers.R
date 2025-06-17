#' @export
connect_db <- function(db_path) {
  DBI::dbConnect(RSQLite::SQLite(), db_path)
}

#' @export
list_tables <- function(con) {
  DBI::dbListTables(con)
}