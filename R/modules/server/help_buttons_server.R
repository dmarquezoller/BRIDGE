source("R/global.R")


help_buttons <- function(input, output, session) {

    general_help <- paste(readLines("./data/general_help.txt"), collapse = "\n")

    observeEvent(input$general_help, {
        showModal(modalDialog(
            title = "General Information",
            HTML(general_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    data_help <- paste(readLines("./data/data_help.txt"), collapse = "\n")

    observeEvent(input$data_help, {
        showModal(modalDialog(
            title = "Data Import",
            HTML(data_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    individual_help <- paste(readLines("./data/individual_help.txt"), collapse = "\n")

    observeEvent(input$individual_help, {
        showModal(modalDialog(
            title = "Individual Exploration",
            HTML(individual_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    raw_int_help <- paste(readLines("./data/raw_int_help.txt"), collapse = "\n")

    observeEvent(input$raw_int_help, {
        showModal(modalDialog(
            title = "Raw Data Integration",
            HTML(raw_int_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    processed_int_help <- paste(readLines("./data/processed_int_help.txt"), collapse = "\n")

    observeEvent(input$processed_int_help, {
        showModal(modalDialog(
            title = "Processed Data Integration",
            HTML(processed_int_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })
}