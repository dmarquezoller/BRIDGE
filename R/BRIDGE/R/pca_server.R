#' @export
pcaServer <- function(id, rv, tbl_name) {
  moduleServer(id, function(input, output, session) {

    # Optional: background task (keeps UI snappy if DESeqTransform is heavy)
    pca_task <- ExtendedTask$new(function(dep) {
      promises::future_promise({
        DESeq2::DESeqTransform(dep)
      })
    })

    observeEvent(input$compute, {
      req(rv$dep_output[[tbl_name]])
      pca_task$invoke(isolate(rv$dep_output[[tbl_name]]))
    }, ignoreInit = TRUE)

    output$plot <- renderPlot({
      res <- pca_task$result()
      req(res)
      DESeq2::plotPCA(res) +
        ggplot2::ggtitle(paste("PCA for", tbl_name)) +
        ggplot2::theme_minimal()
    })
  })
}