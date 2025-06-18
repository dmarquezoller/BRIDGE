#' @export
help_buttons <- function(input, output, session) {

    general_help <- paste(readLines(system.file("helptext", "general_help.txt", package="BRIDGE")), collapse = "\n")

    shiny::observeEvent(input$general_help, {
        showModal(modalDialog(
            title = "General Information",
            shiny::HTML(general_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    data_help <- paste(readLines(system.file("helptext", "data_help.txt", package="BRIDGE")), collapse = "\n")

    shiny::observeEvent(input$data_help, {
        showModal(modalDialog(
            title = "Data Import",
            shiny::HTML(data_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    individual_help <- paste(readLines(system.file("helptext", "individual_help.txt", package="BRIDGE")), collapse = "\n")

    shiny::observeEvent(input$individual_help, {
        showModal(modalDialog(
            title = "Individual Exploration",
            shiny::HTML(individual_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    raw_int_help <- paste(readLines(system.file("helptext", "raw_int_help.txt", package="BRIDGE")), collapse = "\n")

    shiny::observeEvent(input$raw_int_help, {
        showModal(modalDialog(
            title = "Raw Data Integration",
            shiny::HTML(raw_int_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })

    processed_int_help <- paste(readLines(system.file("helptext", "processed_int_help.txt", package="BRIDGE")), collapse = "\n")

    shiny::observeEvent(input$processed_int_help, {
        showModal(modalDialog(
            title = "Processed Data Integration",
            shiny::HTML(processed_int_help), 
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    })
}