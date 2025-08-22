#' @export
VolcanoServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        volcano_ready <- reactiveVal(FALSE)
        last_params <- reactiveVal(NULL) # <- only changes on button click

        volcano_task <- ExtendedTask$new(function(deps) {
            promises::future_promise(
                {
                    dep_output <- deps$dep_output
                    contrast <- deps$contrast
                    datatype <- deps$datatype
                    pcut <- deps$pcut
                    fccut <- deps$fccut
                    highlight <- deps$highlight

                    lfc_col <- paste0(contrast, "_diff")
                    pval_col <- paste0(contrast, "_p.adj")
                    # sig table
                    sig_pg <- DEP2::add_rejections(dep_output, alpha = pcut, lfc = fccut)
                    sig_se <- DEP2::get_signicant(sig_pg, contrast)
                    df_table <- DEP2::get_results(sig_se)

                    # volcano plot df
                    rd <- SummarizedExperiment::rowData(dep_output)

                    # names in the SAME space as the plot
                    name <- switch(datatype,
                        proteomics        = stringr::str_to_title(rd$Gene_Name),
                        phosphoproteomics = rd$pepG,
                        rnaseq            = rownames(rd)
                    )

                    df <- data.frame(
                        name = as.character(name),
                        log2FC = as.numeric(rd[[lfc_col]]),
                        pval = as.numeric(rd[[pval_col]]),
                        stringsAsFactors = FALSE
                    )

                    # normalize highlight into df$name space
                    if (length(highlight)) {
                        if (any(grepl("[,;\\s]", highlight))) {
                            highlight <- unique(unlist(strsplit(highlight, "[,;\\s]+")))
                        }
                        if (identical(datatype, "proteomics")) {
                            highlight <- stringr::str_to_title(highlight)
                        }
                        # if user typed with different case, match case-insensitively:
                        if (!all(highlight %in% df$name)) {
                            map <- match(tolower(highlight), tolower(df$name))
                            highlight <- unique(na.omit(df$name[map]))
                        } else {
                            highlight <- unique(highlight)
                        }
                    } else {
                        highlight <- character(0)
                    }

                    # pre-compute helpers
                    df$neglog10p <- -log10(df$pval)
                    is_sig <- is.finite(df$log2FC) & is.finite(df$neglog10p) &
                        df$pval <= pcut & abs(df$log2FC) >= fccut
                    df$dir <- ifelse(is_sig & df$log2FC > 0, "up",
                        ifelse(is_sig & df$log2FC < 0, "down", "ns")
                    )

                    # KEEP all sig + ALL highlighted, then thin remaining ns
                    if (nrow(df) > 10000) {
                        keep_always <- (df$dir != "ns") | (df$name %in% highlight)
                        ns_idx <- which(!keep_always) # only truly droppable ns
                        if (length(ns_idx)) {
                            ns_keep <- ns_idx[seq(1, length(ns_idx), by = 3L)]
                            keep <- keep_always
                            keep[ns_keep] <- TRUE
                            df <- as.data.frame(df[keep, , drop = FALSE])
                        }
                    }

                    list(
                        df = df, table = df_table, pcut = pcut, fccut = fccut,
                        datatype = datatype, highlight = highlight
                    )
                },
                seed = TRUE
            )
        })

        # Populate highlight choices
        observe({
            req(rv$tables[[tbl_name]], rv$datatype[[tbl_name]])
            choices <- switch(rv$datatype[[tbl_name]],
                phosphoproteomics = rv$tables[[tbl_name]]$pepG,
                rnaseq = if (!is.null(rv$tables[[tbl_name]]$Gene_ID)) {
                    rv$tables[[tbl_name]]$Gene_ID
                } else {
                    rownames(rv$tables[[tbl_name]])
                },
                rv$tables[[tbl_name]]$Gene_Name # proteomics default
            )
            updateSelectizeInput(session, "volcano_search", choices = sort(unique(choices)), server = TRUE)
        })

        # CLICK to compute & freeze params (only-on-click rendering)
        observeEvent(input$compute_volcano,
            {
                req(rv$dep_output[[tbl_name]])
                volcano_ready(TRUE)

                params <- list(
                    contrast   = isolate(input$comparison_volcano),
                    pcut       = isolate(input$volcano_pcutoff), # raw FDR (0..1)
                    fccut      = isolate(input$volcano_fccutoff),
                    highlight  = isolate(input$volcano_search) |> stringr::str_trim(),
                    dep_output = isolate(rv$dep_output[[tbl_name]]),
                    datatype   = isolate(rv$datatype[[tbl_name]])
                )
                last_params(params)
                volcano_task$invoke(params)
            },
            ignoreInit = TRUE
        )

        output$volcano_slot <- renderUI({
            if (!volcano_ready()) {
                return(div(style = "padding:10px; color:#777;", "Click “Compute Volcano” to generate the plot."))
            }
            shinycssloaders::withSpinner(
                plotly::plotlyOutput(session$ns("volcano"), height = "520px"),
                type = 8, color = "#2b8cbe", caption = "Loading..."
            )
        })

        output$volcano <- plotly::renderPlotly({
            req(last_params())
            res <- volcano_task$result()
            req(res)

            df <- res$df

            df <- df[is.finite(df$log2FC) & is.finite(df$pval), , drop = FALSE]
            if (!nrow(df)) {
                showNotification("No finite points to plot for this contrast.", type = "warning")
                return(plotly::plot_ly())
            }
            # message("highlight: ", paste(res$highlight, collapse = ", "))
            # Use lab = df$name for all, and selectLab = highlights to force labels

            p <- EnhancedVolcano::EnhancedVolcano(
                df,
                lab = df$name,
                selectLab = res$highlight,
                x = "log2FC",
                y = "pval",
                title = NULL, subtitle = NULL, caption = NULL,
                pCutoff = res$pcut,
                FCcutoff = res$fccut,
                xlab = "Log[2] fold change",
                ylab = "-Log[10] Padj",
                drawConnectors = FALSE,
                pointSize = 1.8,
                labSize = 3,
                legendPosition = "none"
            ) + ggplot2::aes(text = name)

            # message("Rendering volcano plot with ", nrow(df), " points.", str(p))

            plotly::ggplotly(p + aes(key = name, x = log2FC, y = -log10(pval), tooltip = "text", dynamicTicks = FALSE)) #|>

            # plt <- tryCatch({
            #  plotly::ggplotly(p, tooltip = "text", dynamicTicks = FALSE) #|>
            #  # plotly::toWebGL() |>
            #  # plotly::partial_bundle()
            # }, error = function(e) {
            #  showNotification(
            #    paste("Interactive rendering failed; falling back to basic scatter:", e$message),
            #    type = "error"
            #  )
            #  NULL
            # })
            # if (is.null(plt)) {
            #  return(
            #    plotly::plot_ly(
            #      x = df$log2FC, y = -log10(df$pval),
            #      type = "scatter", mode = "markers",
            #      text = df$name,
            #      hovertemplate = "<b>%{text}</b><br>log2FC=%{x:.2f}<br>-log10(p)=%{y:.2f}<extra></extra>"
            #    )
            #  )
            # }
            # plt
        })

        output$volcano_sig_table <- DT::renderDT({
            res <- volcano_task$result()
            req(res)
            DT::datatable(res$table,
                extensions = "Buttons",
                options = list(
                    scrollX = TRUE, pageLength = 10,
                    dom = "Bfrtip", buttons = c("copy", "csv", "excel", "pdf", "print")
                )
            )
        })
    })
}
