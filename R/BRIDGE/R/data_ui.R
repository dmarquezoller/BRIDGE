#' @export
data_ui <- function(id) {
  ns <- NS(id)
  shiny::fluidRow(
    shiny::uiOutput(ns("all_tables_ui"))
  )
}

#' @export
data_ui_object <- shinydashboard::box(
        title = "Data Selection", width = 12, status = "success", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
        shiny::tags$div(style = "height:8px;"),
       shiny::div(style = "text-align: center;",
            shinyWidgets::actionBttn("data_help", "Help", color = "success" ,icon=shiny::icon("question-circle"), size="sm", style = "bordered")
        ),
        shiny::tags$div(style = "height:8px;"),
        shiny::helpText("Select Species, Table and the desired columns to load the data "),
        shiny::tags$div(style = "height:10px;"),
        shiny::selectInput("species", "Select Species", choices = c("Zebrafish"), multiple = T),
        
        shiny::helpText("At the moment the accepted species are: Zebrafish"),
        shiny::tags$div(style = "height:10px;"),
        shiny::uiOutput("table_selector"),
        shiny::tags$div(style = "height:10px;"),
        shiny::uiOutput("column_selector"),
        shiny::tags$div(style = "height:10px;"),
       shiny::div(style = "text-align: center;",
            shinyWidgets::actionBttn("load_data", shiny::span("Load Data", style = "color: black;"), icon = shiny::span(shiny::icon("play"), style = "color: black;"), style = "jelly", size="sm", class = "btn-success"),
            shiny::br(), shiny::br(),  
            shiny::tags$div(
              style = "text-align: left; padding-left: 10px;",
              shiny::tags$strong("Loaded datasets:"),
              shiny::uiOutput("loaded_info")
            )
        ),
        shiny::tags$div(style = "height:10px;"),
        shiny::selectInput("remove", "Select Table to remove", choices = NULL, selected = NULL, multiple = TRUE),
        shiny::helpText("Write the name of an already loaded table to delete it"),
        shiny::tags$div(style = "height:10px;"),
       shiny::div(style = "text-align: center;",
            shinyWidgets::actionBttn("delete_data", shiny::span("Remove Data", style = "color: black;"), icon = shiny::span(shiny::icon("trash"), style = "color: black;"), style = "jelly", size="sm", class = "btn-success"),
        )       
      )

