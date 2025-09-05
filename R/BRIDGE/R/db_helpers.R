#' @export
connect_db <- function(db_path) {
    DBI::dbConnect(RSQLite::SQLite(), db_path)
}

# Make a pool (for SQLite example)
connect_pool <- function(db_path) {
    pool::dbPool(RSQLite::SQLite(), dbname = db_path)
}

#' @export
list_tables <- function(con) {
    DBI::dbListTables(con)
}
