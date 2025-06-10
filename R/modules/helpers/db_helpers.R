library(DBI)
library(RSQLite)

connect_db <- function(db_path) {
  dbConnect(SQLite(), db_path)
}

list_tables <- function(con) {
  dbListTables(con)
}