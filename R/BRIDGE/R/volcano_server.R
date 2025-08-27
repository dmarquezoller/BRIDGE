#' @export
VolcanoServer <- function(id, rv, cache, tbl_name) {
    moduleServer(id, function(input, output, session) {
        volcano_ready <- reactiveVal(FALSE)
        last_params <- reactiveVal(NULL) # <- only changes on button click
        depflt_cache <- reactiveValues() # in-memory cache of filtered object per (table, columns_key, p, lfc)

        strip_sig <- function(dep_obj) {
            # message("Stripping old sig cols")
            if (methods::is(dep_obj, "DEGdata")) {
                # message("DEGdata: ", str(dep_obj))
                tr <- dep_obj@test_result
                tr <- tr[, !grepl("_significant$|^significant$", colnames(tr))]
                dep_obj@test_result <- tr
            } else {
                # message("SummarizedExperiment: ", class(dep_obj))
                rd <- SummarizedExperiment::rowData(dep_obj)
                rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
                # message("RD: ", class(rd))
                SummarizedExperiment::rowData(dep_obj) <- rd
            }
            # message("Stripped old sig cols")
            dep_obj
        }

        get_dep_result <- function() {
            key <- rv$current_dep_volcano_key[[tbl_name]]
            #message("Prepping DEP results for key: ", key, ", is cached?: ", cache$exists(key))
            if (!is.null(key) && cache$exists(key)) {
                #message("Getting CACHED DEP results for key: ", key)
                return(cache$get(key))
            }
            volcano_task$result()
        }

        get_depflt <- function(params) {
            key <- paste(
                tbl_name, params$columns_key,
                sprintf("pcut=%.4f", params$p_cut),
                sprintf("lfc=%.3f", params$lfc_cut),
                "depflt",
                sep = "_"
            )            
            if (!is.null(depflt_cache[[key]])) {
                return(depflt_cache[[key]])
            }
            # message("Getting DEPflt results for key: ", key)
            dep_output <- params$dep_output
            dep_output <- strip_sig(dep_output)
            dep_flt <- DEP2::add_rejections(dep_output, alpha = params$p_cut, lfc = params$lfc_cut)
            depflt_cache[[key]] <- dep_flt
            dep_flt
        }

        volcano_task <- ExtendedTask$new(function(args) {
            promises::future_promise({                
                dep_output <- args$dep_output
                contrast   <- args$contrast
                datatype   <- args$datatype
                p_cut      <- args$p_cut
                lfc_cut    <- args$lfc_cut
                highlight  <- args$highlight
                #message("DEP RAW: ", class(dep_output))
                lfc_col <- paste0(contrast, "_diff")
                pval_col <- paste0(contrast, "_p.adj")
                # message("DEP: ", class(dep_output))
                dep_output <- strip_sig(dep_output)
                # message("Stripped: ", class(dep_output))
                dep_flt <- DEP2::add_rejections(dep_output, alpha = p_cut, lfc = lfc_cut)
                # message("Filtered: ", class(dep_flt))
                # message("Gene Info: ", head(gene_info), "\nROWS: ", nrow(gene_info), "\nColumns: ", colnames(gene_info))
                # message("Columns: ", paste(colnames(rd), collapse = ", "))
                df <- as.data.frame(SummarizedExperiment::rowData(dep_flt))
                df_table <- df
                message("DF_table: ", head(df_table, 3), "\ndpflt: ", str(SummarizedExperiment::colData(dep_flt)))

                names <- switch(datatype,
                    proteomics        = paste0(stringr::str_to_title(df$Gene_Name), "_", df$Gene_ID),
                    phosphoproteomics = paste0(df$Gene_Name, "_", df$pepG),
                    rnaseq            = rownames(df)
                )
                message("NAMING: ", head(names))

                df <- data.frame(
                    name = names,
                    log2FC = as.numeric(df[[lfc_col]]),
                    pval = as.numeric(df[[pval_col]]),
                    stringsAsFactors = FALSE
                )                
                df <- stats::na.omit(df)
                message("Volcano DF: ", colnames(df), " | ", paste(head(df$name), collapse = ", "))

                list(
                    df = df, table = df_table, pcut = p_cut, fccut = lfc_cut,
                    datatype = datatype, highlight = highlight
                )
            }, seed = TRUE)
        })

        # Populate highlight choices
        observe({
            req(rv$tables[[tbl_name]], rv$datatype[[tbl_name]])
            choices <- switch(
                rv$datatype[[tbl_name]],
                phosphoproteomics = paste0(rv$tables[[tbl_name]]$Gene_Name, "_", rv$tables[[tbl_name]]$pepG),
                rnaseq = if (!is.null(rv$tables[[tbl_name]]$Gene_ID)) {
                            paste0(rv$tables[[tbl_name]]$Gene_Name, "_", rv$tables[[tbl_name]]$Gene_ID)
                         } else {
                            rownames(rv$tables[[tbl_name]])
                         },
                         paste0(rv$tables[[tbl_name]]$Gene_Name, "_", rv$tables[[tbl_name]]$Gene_ID) # proteomics default
            )
            updateSelectizeInput(session, "volcano_search", choices = sort(unique(choices)), server = TRUE)
        })

        # CLICK to compute & freeze params (only-on-click rendering)
        observeEvent(input$compute_volcano,
            {
                req(rv$dep_output[[tbl_name]])
                volcano_ready(TRUE)

                params <- list(
                    contrast     = isolate(input$comparison_volcano),
                    p_cut        = isolate(input$volcano_pcutoff), # raw FDR (0..1)
                    lfc_cut      = isolate(input$volcano_fccutoff),
                    highlight    = isolate(input$volcano_search) |> stringr::str_trim(),
                    dep_output   = isolate(rv$dep_output[[tbl_name]]),
                    datatype     = isolate(rv$datatype[[tbl_name]]),
                    columns_key  = paste(sort(trimws(isolate(rv$time_cols[[tbl_name]]))), collapse = "_")

                )
                last_params(params)
                key <- paste(
                    tbl_name, params$columns_key,
                    sprintf("pcut=%.4f", params$p_cut),
                    sprintf("lfc=%.3f", params$lfc_cut),
                    "dep_volcano_task",
                    sep = "_"
                )
                rv$current_dep_volcano_key[[tbl_name]] <- key

                if (!cache$exists(key)) {
                    volcano_task$invoke(list(
                        dep_output   = params$dep_output,
                        p_cut        = params$p_cut,
                        lfc_cut      = params$lfc_cut,
                        highlight    = params$highlight,
                        contrast     = params$contrast,
                        datatype     = params$datatype
                    ))
                }
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
            params <- req(last_params())
            res <- get_dep_result()
            req(res)

            df <- res$df
            # message("DF: ", colnames(df), " | ", paste(head(df$name), collapse = ", "))
            highlight <- params$highlight
            p_cut <- params$p_cut
            lfc_cut <- params$lfc_cut
            # normalize highlight into df$name space
            # message("Trying to highlight: ", paste(highlight, collapse = ", "), "in ", paste(head(df$name), collapse = ", "))
            if (length(highlight)) {
                if (!all(highlight %in% df$name)) { # if user typed with different case, match case-insensitively
                    map <- match(tolower(highlight), tolower(df$name))
                    highlight <- unique(na.omit(df$name[map]))
                } else {
                    highlight <- unique(highlight)
                }
            } else {
                highlight <- character(0)
            }
            message("Highlighting: ", paste(highlight, collapse = ", "), "\nDFnames: ", paste(head(df$name), collapse = ", "))
            # pre-compute helpers
            df$neglog10p <- -log10(df$pval)
            is_sig <- is.finite(df$log2FC) & is.finite(df$neglog10p) &
                df$pval <= p_cut & abs(df$log2FC) >= lfc_cut
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
            df <- df[is.finite(df$log2FC) & is.finite(df$pval), , drop = FALSE]
            if (!nrow(df)) {
                showNotification("No finite points to plot for this contrast.", type = "warning")
                return(plotly::plot_ly())
            }
            message("highlight: ", paste(highlight, collapse = ", "))
            # Use lab = df$name for all, and selectLab = highlights to force labels

            keyvals <- ifelse(
                tolower(df$name) %in% tolower(highlight), "black",
                ifelse(df$log2FC < -1 * lfc_cut, "royalblue",
                    ifelse(df$log2FC > lfc_cut, "red",
                        "grey80"
                    )
                )
            )
            keyvals[is.na(keyvals)] <- "grey80"
            names(keyvals)[keyvals == "black"] <- "highlight"
            names(keyvals)[keyvals == "grey80"] <- "mid"
            names(keyvals)[keyvals == "royalblue"] <- "low"
            names(keyvals)[keyvals == "red"] <- "high"

            # message(length(keyvals), " points; highlights: ", length(res$highlight), " (", paste(res$highlight, collapse = ", "), ")")

            p <- EnhancedVolcano::EnhancedVolcano(
                df,
                lab = df$name,
                selectLab = highlight,
                x = "log2FC",
                y = "pval",
                title = NULL, subtitle = NULL, caption = NULL,
                pCutoff = p_cut,
                FCcutoff = lfc_cut,
                colCustom = keyvals, 
                xlab = "Log[2] fold change",
                ylab = "-Log[10] Padj",
                drawConnectors = FALSE,
                pointSize = 1.8,
                labSize = 3,
                legendPosition = "none"
            ) + ggplot2::aes(text = name)

            # message("Rendering volcano plot with ", nrow(df), " points.", str(p))

            plotly::ggplotly(p + aes(key = name, x = log2FC, y = -log10(pval), tooltip = "text", dynamicTicks = FALSE)) #|>
        })

        output$volcano_sig_table <- DT::renderDT({
           params <- req(last_params())
           res <- get_dep_result()
           req(res)
           dep_flt <- get_depflt(params)
           # message("Filtered DEGs: ", class(dep_flt))

           key <- rv$current_dep_volcano_key[[tbl_name]]
           if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

           dep_flt <- get_depflt(params)
           if (methods::is(dep_flt, "DEGdata")) {
               df <- res$table
               df$name <- rownames(df)
               df$Gene_ID <- gsub("^(.*)_", "", rownames(df))
               gene_map <- rv$tables[[tbl_name]][, c("Gene_ID", "Gene_Name")]
               df <- dplyr::left_join(df, gene_map, by = "Gene_ID")

               sig <- as.data.frame(dep_flt@test_result)
               sig$Gene_ID <- rownames(dep_flt@test_result)
               sig <- dplyr::left_join(sig, gene_map, by = "Gene_ID")
               sig_genes <- sig$Gene_Name[sig$significant]

               df_filtered <- df[df$Gene_Name %in% sig_genes, , drop = FALSE]
               rownames(df_filtered) <- NULL
               df_filtered <- df_filtered %>% select(Gene_Name, Gene_ID, name, everything())
           } else {
               rd <- SummarizedExperiment::rowData(dep_flt)
               sig_genes <- rd$Gene_Name[rd$significant]
               df_filtered <- res$table[res$table$Gene_Name %in% sig_genes, , drop = FALSE] %>% dplyr::select(-ID)
               df_filtered$name <- rownames(df_filtered)
               rownames(df_filtered) <- NULL
           }

           DT::datatable(df_filtered,
               extensions = "Buttons",
               options = list(
                   scrollX = TRUE, pageLength = 10,
                   dom = "Bfrtip", buttons = c("copy", "csv", "excel", "pdf", "print")
               )
           )
        })
    })
}
