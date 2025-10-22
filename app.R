if (!require(BRIDGE)) {
    library(devtools)
    devtools::document("R/BRIDGE")
    devtools::install("R/BRIDGE", keep_source = T, upgrade = "never")    
}
print("Loading BRIDGE package")
library(BRIDGE)
print("Loaded BRIDGE package")

suppressPackageStartupMessages({
    library(future.callr) # or future for multisession
    library(promises)
    library(tidyverse)
    library(shiny)
    library(shinybusy)
    library(DEP2)
    library(NbClust)
    library(SummarizedExperiment)
    library(ggplot2)
    library(plotly)
    library(ComplexHeatmap)
    library(matrixStats)
    library(gtools)
})
# library(pool)
# library(shinydashboard)
# library(shinyWidgets)
ht_opt$message <- FALSE
# future::plan(multisession, workers = 4)
future::plan(future.callr::callr, workers = 4)
set.seed(42)

print("Getting CLI arguments")

# Determine database path
get_db_path <- function() {
    # 1. Try environment variable (for Docker/Shiny Server)
    db_path <- Sys.getenv("BRIDGE_DB_PATH", unset = NA)
    if (!is.na(db_path) && file.exists(db_path)) return(db_path)
    # 2. Try default Docker path
    if (file.exists("/srv/data/database.db")) return("/srv/data/database.db")
    # 3. If interactive, prompt user
    if (interactive()) return(file.choose())
    # 4. If command line, use first argument
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0 && file.exists(args[1])) return(args[1])
    stop("No valid database path found.")
}
db_path <- get_db_path()
print(paste("db_path:", db_path))

get_port <- function() {
    # 1. Try environment variable (for Docker/Shiny Server)
    port_env <- Sys.getenv("BRIDGE_PORT", unset = NA)
    if (!is.na(port_env) && nzchar(port_env)) return(as.integer(port_env))
    # 2. Try command line argument
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 1 && !is.na(as.integer(args[2]))) return(as.integer(args[2]))
    # 3. Default
    return(3838L)
}
port <- get_port()
print(paste("port:", port))


print("About to define ui and server")

# Define UI and server
ui <- BRIDGE::ui
print("UI")
print(exists("ui"))
#print(exists("ui", where=asNamespace("BRIDGE")))

server <- function(input, output, session) {
    BRIDGE::server_function(input, output, session, db_path)
}
print("SERVER")
print(exists("server"))

shiny::runApp(
    shiny::shinyApp(ui = ui, server = server), port = port
)
# Only run the app if not under Shiny Server
#if (interactive() || (Sys.getenv("RUN_STANDALONE", "0") == "1")) {
#    shiny::runApp(
#        shiny::shinyApp(ui = ui, server = server),
#        port = port
#    )
#}