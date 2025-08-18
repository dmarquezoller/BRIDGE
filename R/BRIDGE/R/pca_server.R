#' @export
pcaServer <- function(id, rv) 
    moduleServer(
    id,
    function(input, output, session){
    shiny::observe({
        lapply(rv$table_names, function(tbl_name) {
            shiny::observeEvent(input[[paste0("compute_pca_", tbl_name)]], {
                output[[paste0("pca_", tbl_name)]] <- shiny::renderPlot({
                    dep_output <- rv$dep_output[[tbl_name]]
                    pca_plot_counts <- DESeq2::DESeqTransform(dep_output)
                    DESeq2::plotPCA(pca_plot_counts)
                })
            })
        })
    })
}
)