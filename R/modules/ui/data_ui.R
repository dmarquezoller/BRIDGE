data_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    uiOutput(ns("all_tables_ui"))
  )
}

data_ui_object <- box(
        title = "Data Selection", width = 12, status = "success", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
        tags$div(style = "height:8px;"),
        div(style = "text-align: center;",
            actionBttn("data_help", "Help", color = "success" ,icon=icon("question-circle"), size="sm", style = "bordered")
        ),
        tags$div(style = "height:8px;"),
        helpText("Select Species, Table and the desired columns to load the data "),
        tags$div(style = "height:10px;"),
        selectInput("species", "Select Species", choices = c("Zebrafish"), multiple = T),
        
        helpText("At the moment the accepted species are: Zebrafish"),
        tags$div(style = "height:10px;"),
        uiOutput("table_selector"),
        tags$div(style = "height:10px;"),
        uiOutput("column_selector"),
        tags$div(style = "height:10px;"),
        div(style = "text-align: center;",
            actionBttn("load_data", span("Load Data", style = "color: black;"), icon = span(icon("play"), style = "color: black;"), style = "jelly", size="sm", class = "btn-success"),
            br(), br(),  
            tags$div(
              style = "text-align: left; padding-left: 10px;",
              tags$strong("Loaded datasets:"),
              uiOutput("loaded_info")
            )
        ),
        tags$div(style = "height:10px;"),
        selectInput("remove", "Select Table to remove", choices = NULL, selected = NULL, multiple = TRUE),
        helpText("Write the name of an already loaded table to delete it"),
        tags$div(style = "height:10px;"),
        div(style = "text-align: center;",
            actionBttn("delete_data", span("Remove Data", style = "color: black;"), icon = span(icon("trash"), style = "color: black;"), style = "jelly", size="sm", class = "btn-success"),
        )       
      )

