library(BRIDGE)
library(future)
library(shiny)
library(shinydashboard)
library(shinyWidgets)

pdf(file = NULL)

future::plan(multisession)

# Function to run the app with a command-line argument for db_path
bridge <- function() {
  # Retrieve command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  # Check if db_path is provided
  if (length(args) < 1) {
    stop("Error: Please provide the database path as a command-line argument.")
  }
  
  db_path <- args[1]# First argument is the database path
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
}

# Run the app
bridge()
