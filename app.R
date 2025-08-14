if(!require(BRIDGE)){
  library(devtools)
  devtools::document("R/BRIDGE")
  devtools::install("R/BRIDGE", keep_source = T, upgrade = "never")
  library(BRIDGE)  
}
#library(future)
library(future.callr)  # or future for multisession
library(promises)
library(tidyverse)
library(shiny)
#library(shinydashboard)
#library(shinyWidgets)
options(shiny.fullstacktrace = TRUE, shiny.sanitize.errors = FALSE)

shiny_error_trap <- function(expr) {
  withCallingHandlers(expr, error = function(e) {
    cat("\n=== HARD TRACE ===\n")
    cs <- sys.calls()
    for (i in seq_along(cs)) cat(sprintf("%2d: %s\n", i, deparse(cs[[i]], nlines = 1L)))
    cat("=== /HARD TRACE ===\n\n")
    stop(e)
  })
}

trace(jsonlite:::asJSON, tracer = quote({
  if (is.environment(x)) {
    cat("\n>>> asJSON got an ENVIRONMENT. Call stack:\n")
    print(sys.calls())
    cat("<<<\n")
  }
}), print = FALSE)

# pdf(file = NULL)  # delete, can trigger Warning: CallrFuture (<unnamed-6>) added, removed, or modified devices. A future expression must close any opened devices and must not close devices it did not open. Details: 1 devices differ: index=2, before=‘NA’, after=‘pdf’. See also help("future.options", package = "future")
# # 1 thread per worker to allow true parallelism
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
library(future); library(future.callr)
options(
  #future.debug = TRUE,  # Enable debugging for future  
  future.seed = TRUE
)
future::plan(future.callr::callr, workers = 2L)
#message("workers = ", future::nbrOfWorkers())

set.seed(42)
# Function to run the app with a command-line argument for db_path
bridge <- function() {
  shiny_error_trap({
    # Retrieve command-line arguments
    args <- commandArgs(trailingOnly = TRUE)
    
    # Handle database path argument for both terminal and interactive use
    if (interactive()) {
      if (length(args) == 0) {
        # Prompt user to select the database file interactively
        db_path <- file.choose()
      } else {
        db_path <- args[1]
      }
    } else {
      if (length(args) == 0) {
        stop("Please provide the path to the database as an argument.")
      }
      db_path <- args[1]
    }
    
    # Run the Shiny app, passing db_path to the server function
    app <- shiny::shinyApp(
      ui = BRIDGE::ui, #ui retrieved from initial_ui.R
      server = function(input, output, session) {
        server_function(input, output, session, db_path) #Function in general_server.R
      }
    )
    if (length(args) > 1 ) {
      port <- as.integer(args[2]) # Second argument is the port number
      if (is.na(port) || port <= 0 || port > 65535) {
        stop("Error: Invalid port number. Please provide a valid port number between 1 and 65535.")
      }
      shiny::runApp(app, port = port) # Run the app on the specified port
    }
    else {
      # If no port is specified, run on the default port (3838)
      shiny::runApp(app, port = 3838)
    }
  })
}

# Run the app
bridge()
