#' @export
ui <- shinydashboard::dashboardPage(
    skin = "blue",
    shinydashboard::dashboardHeader(title = "BRIDGE"),
    shinydashboard::dashboardSidebar(
        shinydashboard::sidebarMenu(
            shinydashboard::menuItem("Individual Exploration", tabName = "raw_data", icon = shiny::icon("table")),
            shinydashboard::menuItem("Integration", tabName = "Integration", icon = shiny::icon("th")),
            shinyWidgets::actionBttn("general_help", "Help", icon = shiny::icon("question-circle"), size = "sm", style = "bordered")
        )
    ),
    shinydashboard::dashboardBody(
        # Busy message overlay
        shinybusy::add_busy_spinner(spin = "radar", position = "bottom-right", timeout = 200),
        # Custom font styling for the dashboard title
        shiny::tags$head(
            shiny::tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css?family=Inknut+Antiqua"),
            shiny::tags$style(shiny::HTML("
        .main-header .logo {
          font-family: 'Inknut Antiqua', sans-serif;
        }
      "))
        ),
        shiny::fluidRow(
            shiny::column(
                width = 9,
                shinydashboard::tabItems(
                    shinydashboard::tabItem(
                        tabName = "raw_data",
                        shiny::fluidRow(
                            shiny::uiOutput("all_tables_ui")
                        )
                    ),
                    shinydashboard::tabItem(
                        tabName = "Integration",
                        shiny::fluidRow(
                            shiny::uiOutput("integration_ui")
                        )
                    )
                )
            ),
            shiny::column(
                width = 3,
                data_ui_object
            )
        )
    )
)
