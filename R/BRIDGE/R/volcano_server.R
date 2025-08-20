#' @export
VolcanoServer <- function(id, rv, tbl_name) {
  moduleServer(id, function(input, output, session) {

    volcano_ready <- reactiveVal(FALSE)
    last_params   <- reactiveVal(NULL)   # <- only changes on button click

    volcano_task <- ExtendedTask$new(function(deps) {
      promises::future_promise({
        dep_output <- deps$dep_output
        contrast   <- deps$contrast
        datatype   <- deps$datatype
        pcut       <- deps$pcut      # raw 0..1
        fccut      <- deps$fccut
        highlight  <- deps$highlight

        lfc_col  <- paste0(contrast, "_diff")
        pval_col <- paste0(contrast, "_p.adj")
        rd       <- SummarizedExperiment::rowData(dep_output)

        # 1) compact DF once, with names normalized up-front
        name <- switch(datatype,
          proteomics        = stringr::str_to_title(rd$Gene_Name),
          phosphoproteomics = rd$pepG,
          rnaseq            = rownames(rd)
        )
        df <- data.frame(
          name   = name,
          log2FC = as.numeric(rd[[lfc_col]]),
          pval   = as.numeric(rd[[pval_col]]),
          stringsAsFactors = FALSE
        )
        # pre-compute helpers
        df$neglog10p <- -log10(df$pval)
        is_sig <- is.finite(df$log2FC) & is.finite(df$neglog10p) &
                  df$pval <= pcut & abs(df$log2FC) >= fccut
        df$dir <- ifelse(is_sig & df$log2FC > 0, "up",
                  ifelse(is_sig & df$log2FC < 0, "down", "ns"))

        # 2) significant table (same contrast)
        sig_se   <- DEP2::get_signicant(dep_output, contrast)
        df_table <- DEP2::get_results(sig_se)

        # normalize highlight to df$name space
        if (length(highlight)) {
          if (identical(datatype, "proteomics")) highlight <- stringr::str_to_title(highlight)
          highlight <- intersect(unique(highlight), df$name)
        }

        # 3) optional thinning of massive non-significant cloud (keeps all sig)
        if (nrow(df) > 50000) {
          keep <- which(df$dir != "ns")
          ns   <- which(df$dir == "ns")
          ns_keep <- ns[seq(1, length(ns), by = 3L)]  # keep ~33% ns points
          df <- df[c(keep, ns_keep), , drop = FALSE]
        }

        list(df=df, table=df_table, pcut=pcut, fccut=fccut,
             datatype=datatype, highlight=highlight)
      }, seed = TRUE)
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

    # CLICK to compute & freeze params (only-on-click rendering)
    observeEvent(input$compute_volcano, {
      req(rv$dep_output[[tbl_name]])
      volcano_ready(TRUE)

      params <- list(
        contrast   = isolate(input$comparison_volcano),
        pcut       = isolate(input$volcano_pcutoff),   # raw FDR (0..1)
        fccut      = isolate(input$volcano_fccutoff),
        highlight  = isolate(input$volcano_search) |> stringr::str_trim(),
        dep_output = isolate(rv$dep_output[[tbl_name]]),
        datatype   = isolate(rv$datatype[[tbl_name]])
      )
      last_params(params)
      volcano_task$invoke(params)
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
      req(last_params())                # <- only rerender after the button
      res <- volcano_task$result(); req(res)
      df  <- res$df

      # Only label highlights; everything else unlabeled
      lab_vec <- if (length(res$highlight)) ifelse(df$name %in% res$highlight, df$name, "") else ""

      # Lighter EnhancedVolcano: no connectors, no extra legends
      p <- EnhancedVolcano::EnhancedVolcano(
        df,
        lab         = lab_vec,
        x           = "log2FC",
        y           = "pval",
        title       = NULL,
        subtitle    = NULL,
        caption     = NULL,
        pCutoff     = res$pcut,
        FCcutoff    = res$fccut,
        drawConnectors = FALSE,
        pointSize   = 1.8,     # scalar (faster); highlights will still have labels
        labSize     = 3,
        legendPosition = "none"
      ) + ggplot2::aes(text = name)

      plt <- plotly::ggplotly(p, tooltip = "text", dynamicTicks = FALSE)
      plt <- plotly::toWebGL(plt)          # use WebGL renderer under the hood
      plt <- plotly::partial_bundle(plt)   # ship a lighter JS bundle
      plt
    })

    output$volcano_sig_table <- DT::renderDT({
      res <- volcano_task$result(); req(res)
      DT::datatable(res$table, extensions="Buttons",
                    options=list(scrollX=TRUE, pageLength=10,
                                 dom="Bfrtip", buttons=c('copy','csv','excel','pdf','print')))
    })
  })
}