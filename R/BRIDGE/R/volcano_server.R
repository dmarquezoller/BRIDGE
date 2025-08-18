#' @export
VolcanoServer <- function(id, rv, tbl_name) {
  moduleServer(id, function(input, output, session) {

    volcano_ready <- reactiveVal(FALSE)

    volcano_task <- ExtendedTask$new(function(deps) {
      promises::future_promise({
        dep_output <- deps$dep_output
        contrast   <- deps$contrast
        datatype   <- deps$datatype
        pcut       <- deps$pcut
        fccut      <- deps$fccut
        highlight  <- deps$highlight

        lfc_col  <- paste0(contrast, "_diff")
        pval_col <- paste0(contrast, "_p.adj")
        res      <- SummarizedExperiment::rowData(dep_output)

        if (datatype == "proteomics") {
          df <- data.frame(gene_names = stringr::str_to_title(res$Gene_Name),
                           pval = res[[pval_col]], log2FC = res[[lfc_col]],
                           stringsAsFactors = FALSE)
          sig_se   <- DEP2::get_signicant(dep_output, contrast)
          df_table <- DEP2::get_results(sig_se)
        } else if (datatype == "phosphoproteomics") {
          df <- data.frame(peptide = res$pepG,
                           pval = res[[pval_col]], log2FC = res[[lfc_col]],
                           stringsAsFactors = FALSE)
          sig_se   <- DEP2::get_signicant(dep_output, contrast)
          df_table <- DEP2::get_results(sig_se)
        } else { # rnaseq
          df <- data.frame(gene_ID = rownames(res),
                           pval = res[[pval_col]], log2FC = res[[lfc_col]],
                           stringsAsFactors = FALSE)
          df_table <- DEP2::get_results(DEP2::get_signicant(dep_output, contrast))
        }

        list(df=df, table=df_table, pcut=pcut, fccut=fccut,
             datatype=datatype, highlight=highlight)
      })
    })

    # Populate highlight choices for THIS module's namespace
    observe({
      req(rv$tables[[tbl_name]], rv$datatype[[tbl_name]])
      choices <- if (identical(rv$datatype[[tbl_name]], "phosphoproteomics")) {
        rv$tables[[tbl_name]]$pepG
      } else {
        rv$tables[[tbl_name]]$Gene_Name
      }
      updateSelectizeInput(session, "volcano_search", choices = choices, server = TRUE)
    })

    
    observeEvent(input$compute_volcano, {
      req(rv$dep_output[[tbl_name]])
      volcano_ready(TRUE)

      deps <- list(
        contrast   = isolate(input$comparison_volcano),
        pcut       = isolate(input$volcano_pcutoff),
        fccut      = isolate(input$volcano_fccutoff),
        highlight  = isolate(input$volcano_search) |> stringr::str_trim() |> stringr::str_to_title(),
        dep_output = isolate(rv$dep_output[[tbl_name]]),
        datatype   = isolate(rv$datatype[[tbl_name]])
      )
      volcano_task$invoke(deps)
    }, ignoreInit = TRUE)

    output$volcano_slot <- renderUI({
      if (!volcano_ready()) {
        return(div(style="padding:10px; color:#777;", "Click “Compute Volcano” to generate the plot."))
      }
      shinycssloaders::withSpinner(
        plotly::plotlyOutput(session$ns("volcano"), height = "520px"),
        type = 8, color = "#2b8cbe", caption = "Loading..."
      )
    })

    output$volcano <- plotly::renderPlotly({
      res <- volcano_task$result(); req(res)
      if (res$datatype == "proteomics") {
        text_col <- res$df$gene_names; lab_map <- ggplot2::aes(text = gene_names)
        p <- EnhancedVolcano::EnhancedVolcano(
          res$df, lab=res$df$gene_names, selectLab="a", x="log2FC", y="pval",
          title="", pCutoff=res$pcut, FCcutoff=res$fccut,
          pointSize = ifelse(text_col %in% res$highlight, 3, 1),
          legendPosition="none"
        ) + lab_map + ggplot2::labs(color="Legend")
      } else if (res$datatype == "phosphoproteomics") {
        text_col <- res$df$peptide; lab_map <- ggplot2::aes(text = peptide)
        p <- EnhancedVolcano::EnhancedVolcano(
          res$df, lab=res$df$peptide, selectLab="a", x="log2FC", y="pval",
          title="", pCutoff=res$pcut, FCcutoff=res$fccut,
          pointSize = ifelse(text_col %in% res$highlight, 3, 1),
          legendPosition="none"
        ) + lab_map + ggplot2::labs(color="Legend")
      } else {
        text_col <- res$df$gene_ID; lab_map <- ggplot2::aes(text = gene_ID)
        p <- EnhancedVolcano::EnhancedVolcano(
          res$df, lab=res$df$gene_ID, selectLab="a", x="log2FC", y="pval",
          title="", pCutoff=res$pcut, FCcutoff=res$fccut,
          pointSize = ifelse(text_col %in% res$highlight, 3, 1),
          legendPosition="none"
        ) + lab_map + ggplot2::labs(color="Legend")
      }
      plotly::ggplotly(p + ggplot2::aes(x = log2FC, y = -log10(pval)), tooltip = "text")
    })

    output$volcano_sig_table <- DT::renderDT({
      res <- volcano_task$result(); req(res)
      DT::datatable(res$table, extensions="Buttons",
                    options=list(scrollX=TRUE, pageLength=10,
                                 dom="Bfrtip", buttons=c('copy','csv','excel','pdf','print')))
    })
  })
}