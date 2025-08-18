#' @export
VolcanoUI  <- function(id, tbl_name, contrasts) {   
                        ns <- NS(id)                        
                        shinydashboard::box(title = "Volcano", width = 12, solidHeader = TRUE, status = "info",
                        shiny::fluidRow(),
                        shiny::fluidRow(
                        shinydashboard::box(title = "Settings", width = 5, solidHeader = T, status = "info",
                            shinyWidgets::virtualSelectInput(ns(paste0("comparison_volcano_", tbl_name)), "Select Comparison:", choices = contrasts, selected = contrasts[1]),
                            shiny::selectizeInput(ns(paste0("volcano_search_",tbl_name)), "Search for Gene:", choices = NULL),
                            shinyWidgets::chooseSliderSkin("Flat", color="#3d8dbc"),
                            shiny::numericInput(ns(paste0("volcano_pcutoff_", tbl_name)), "FDR Threshold:", 
                                        min = 0, max = 10, value = 0.05, step = 0.01),
                            shiny::numericInput(ns(paste0("volcano_fccutoff_", tbl_name)), "FC Threshold:", 
                                        min = 0, max = 10, value = 1, step = 0.1),
                            shinyWidgets::actionBttn(ns(paste0("compute_volcano_", tbl_name)), shiny::span("Compute  Volcano", style = "color: white;"), style = "simple", color = "primary", size = "sm")
                        ),
                        shinydashboard::box(
                            title = "Volcano Plot", width = 7, solidHeader = TRUE, status = "info",
                            shinycssloaders::withSpinner(
                            plotly::plotlyOutput(ns(paste0("volcano_", tbl_name))),
                            type = 8,
                            color = "#2b8cbe", 
                            caption = "Loading..."
                            )  
                        )
                        ),
                        h5(),
                        shinydashboard::box(title="Significant Genes", width=12, solidHeader = T, status="info", collapsible = T, collapsed = F,
                            h3(),
                            DT::DTOutput(ns(paste0("volcano_sig_table_", tbl_name))),
                            h3()
                        )
                        )
}
