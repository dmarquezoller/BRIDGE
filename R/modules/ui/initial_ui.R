source("R/modules/ui/data_ui.R")

ui <- dashboardPage(
  skin = "blue", 
  dashboardHeader(title = "BRIDGE"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Raw Data Tables", tabName = "raw_data", icon = icon("table")),
      menuItem("Integration", tabName = "Integration", icon = icon("th")),
      actionBttn("general_help", "Help", icon=icon("question-circle"), size="sm", style = "bordered")
    )
  ),
  dashboardBody(
    # Custom font styling for the dashboard title
    tags$head(
      tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Inknut+Antiqua"),
      tags$style(HTML("
        .main-header .logo {
          font-family: 'Inknut Antiqua', sans-serif;
        }
      "))
    ),
    
    fluidRow(
      column(
        width = 9,
        tabItems(
          tabItem(
            tabName = "raw_data",
            fluidRow(
              uiOutput("all_tables_ui")
            )
          ),
          tabItem(
            tabName = "Integration",
            fluidRow(
              uiOutput("integration_ui")
            )
          )
        )
      ),
      column(
        width = 3,
        data_ui_object
      )
    )
  )
)
