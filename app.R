library(BRIDGE)
library(future)

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
    ui = ui, #ui retrieved from initial_ui.R
    server = function(input, output, session) {
      server_function(input, output, session, db_path) #Function in general_server.R
    }
  )
  shiny::runApp(app)
}

# Run the app
bridge()
