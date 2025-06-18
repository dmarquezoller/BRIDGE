#' @export
pca_server <- function(input, output, session, rv) {
    observe({
        lapply(rv$table_names, function(tbl_name) {
            observeEvent(input[[paste0("compute_pca_", tbl_name)]], {
                output[[paste0("pca_", tbl_name)]] <- renderPlot({
                    dep_output <- rv$dep_output[[tbl_name]]
                    pca_plot_counts <- DESeq2::DESeqTransform(dep_output)
                    plotPCA(pca_plot_counts)
                })
            })
        })
    })
}