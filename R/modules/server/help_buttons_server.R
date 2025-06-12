source("R/global.R")

observeEvent(input$general_help, {
    showModal(modalDialog(
        title = "",
        HTML(), 
        easyClose = TRUE,
        footer = modalButton("Close")
    ))
})