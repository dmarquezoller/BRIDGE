#' @export
heatmap_server <- function(id, rv) {
    moduleServer(id, function(input, output, session) {
        shiny::observe({
            lapply(rv$table_names, function(tbl_name) {
                output[[paste0("raw_ht_", tbl_name)]] <- shiny::renderPlot({
                    ComplexHeatmap::Heatmap(rv$ht_matrix[[tbl_name]], show_row_dend = FALSE, show_row_names = FALSE)
                })
            })
        })
    })
}

#' @export
RawHeatmapServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        ready <- reactiveVal(FALSE)

        task <- ExtendedTask$new(function(raw_input) {
            promises::future_promise({
                raw_data <- raw_input$raw_data
                timepoint_cols <- raw_input$timepoint_cols
                datatype <- raw_input$datatype

                clean <- as.matrix(raw_data[timepoint_cols])
                rn <- switch(datatype,
                    proteomics        = raw_data$Gene_Name,
                    rnaseq            = raw_data$Gene_Name,
                    phosphoproteomics = raw_data$pepG
                )
                rownames(clean) <- rn
                list(clean_matrix = clean)
            })
        })

        observeEvent(input$compute,
            {
                ready(TRUE)
                task$invoke(list(
                    raw_data       = isolate(rv$tables[[tbl_name]]),
                    timepoint_cols = isolate(rv$time_cols[[tbl_name]]),
                    datatype       = isolate(rv$datatype[[tbl_name]])
                ))
            },
            ignoreInit = TRUE
        )

        # Mount spinner only after "Compute" is pressed
        output$plot_slot <- renderUI({
            if (!ready()) {
                return(div(
                    style = "padding:10px; color:#777;",
                    "Click “Compute Raw Heatmap” to generate the plot."
                ))
            }
            shinycssloaders::withSpinner(
                plotOutput(session$ns("plot"), height = "520px"),
                type = 8, color = "#2b8cbe", caption = "Loading..."
            )
        })

        output$plot <- renderPlot({
            res <- task$result()
            req(res)
            ComplexHeatmap::Heatmap(res$clean_matrix,
                show_row_dend = FALSE,
                show_row_names = FALSE
            )
        })
    })
}


#### THIS FUNCTION CONTAINS THE SERVER FOR THE DEP HEATMAP
#' @export
DepHeatmapServer <- function(id, rv, cache, tbl_name) {
    moduleServer(id, function(input, output, session) {
        heatmap_ready <- reactiveVal(FALSE)
        depflt_cache <- reactiveValues() # in-memory cache of filtered object per (table, columns_key, p, lfc)

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

            dep_output <- params$dep_output
            # strip old significant cols once
            if (methods::is(dep_output, "DEGdata")) {
                tr <- dep_output@test_result
                tr <- tr[, !grepl("_significant$|^significant$", colnames(tr))]
                dep_output@test_result <- tr
            } else {
                rd <- SummarizedExperiment::rowData(dep_output)
                rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
                SummarizedExperiment::rowData(dep_output) <- rd
            }

            dep_flt <- DEP2::add_rejections(dep_output, alpha = params$p_cut, lfc = params$lfc_cut)
            depflt_cache[[key]] <- dep_flt
            dep_flt
        }
        last_params <- reactiveVal(NULL) # frozen on click

        # Background task: compute optimal_k from the filtered (significant) matrix + a DF snapshot
        heatmap_task <- ExtendedTask$new(function(args) {
            promises::future_promise(
                {
                    dep_output <- args$dep_output
                    p_cut <- args$p_cut
                    lfc_cut <- args$lfc_cut

                    # filter with current cutoffs
                    dep_flt <- DEP2::add_rejections(dep_output, alpha = p_cut, lfc = lfc_cut)
                    sig <- DEP2::get_signicant(dep_flt) # only significant rows
                    mat <- SummarizedExperiment::assay(sig)
                    mat <- as.matrix(mat)

                    # sanitize matrix (rows/cols need finite variance and ≥2 finite values)
                    row_ok <- rowSums(is.finite(mat)) >= 2 & apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
                    col_ok <- colSums(is.finite(mat)) >= 2 & apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
                    mat <- mat[row_ok, col_ok, drop = FALSE]
                    mat_scaled <- safe_row_scale(mat)

                    # avoid device warnings on workers
                    tmp <- tempfile(fileext = ".pdf")
                    pdf(tmp)
                    on.exit(
                        {
                            dev.off()
                            unlink(tmp)
                        },
                        add = TRUE
                    )

                    optimal_k <- safe_nbclust(mat_scaled, k_min = 2, k_max = 10)
                    if (is.na(optimal_k)) {
                        showNotification("Not enough clean data to estimate k (after filtering).", type = "warning")
                        return(invisible())
                    }

                    # table snapshot based on same filtered/significant set
                    gene_info <- as.data.frame(SummarizedExperiment::rowData(sig))
                    df <- cbind(gene_info, as.data.frame(mat))
                    df <- df[, c(colnames(gene_info), colnames(mat)), drop = FALSE]
                    df <- stats::na.omit(df)

                    list(optimal_k = optimal_k, df = df)
                },
                seed = TRUE
            )
        })

        # Click: freeze params and kick off background compute (key includes cutoffs)
        observeEvent(input$recompute_heatmap,
            {
                req(rv$dep_output[[tbl_name]])
                heatmap_ready(TRUE)

                params <- list(
                    p_cut       = input$heatmap_pcutoff, # raw FDR 0..1
                    lfc_cut     = input$heatmap_fccutoff,
                    clustering  = isTRUE(input$clustering),
                    k           = input$num_clusters,
                    dep_output  = isolate(rv$dep_output[[tbl_name]]),
                    columns_key = paste(isolate(rv$time_cols[[tbl_name]]), collapse = "_")
                )
                last_params(params)

                # heatmap_task only depends on table + selected columns + cutoffs
                key <- paste(
                    tbl_name, params$columns_key,
                    sprintf("pcut=%.4f", params$p_cut),
                    sprintf("lfc=%.3f", params$lfc_cut),
                    "dep_heatmap_task",
                    sep = "_"
                )
                rv$current_dep_heatmap_key[[tbl_name]] <- key

                if (!cache$exists(key)) {
                    heatmap_task$invoke(list(
                        dep_output = params$dep_output,
                        p_cut = params$p_cut, 
                        lfc_cut = params$lfc_cut
                    ))
                }
            },
            ignoreInit = TRUE
        )

        get_dep_result <- function() {
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && cache$exists(key)) {
                return(cache$get(key))
            }
            heatmap_task$result()
        }

        # Spinner container
        output$ht_slot <- renderUI({
            if (!heatmap_ready()) {
                return(div(
                    style = "padding:10px; color:#777;",
                    "Choose settings and click “Recompute heatmap”."
                ))
            }
            shinycssloaders::withSpinner(
                plotOutput(session$ns("ht"), height = "520px"),
                type = 8, color = "#2b8cbe", caption = "Loading..."
            )
        })

        # Plot (depends only on frozen params from the click)
        output$ht <- renderPlot({
            params <- req(last_params())
            res <- get_dep_result()
            req(res)
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

            dep_flt <- get_depflt(params)
            if (params$clustering) {
                k <- as.integer(params$k)
                k_max <- max(2L, nrow(dep_flt) - 1L)
                if (!is.finite(k) || k < 2L) k <- 2L
                if (k > k_max) {
                    showNotification(sprintf(
                        "Reducing k from %d to %d (only %d rows).",
                        input$num_clusters, k_max, nrow(dep_flt)
                    ), type = "warning")
                    k <- k_max
                }
                DEP2::plot_heatmap(dep_flt, kmeans = TRUE, k = k)
            } else {
                DEP2::plot_heatmap(dep_flt)
            }
        })

        output$ht_sig <- DT::renderDT({
            params <- req(last_params())
            res <- get_dep_result()
            req(res)
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

            dep_flt <- get_depflt(params)
            if (methods::is(dep_flt, "DEGdata")) {
                df <- res$df
                df$Gene_ID <- rownames(df)
                gene_map <- rv$tables[[tbl_name]][, c("Gene_ID", "Gene_Name")]
                df <- dplyr::left_join(df, gene_map, by = "Gene_ID")

                sig <- as.data.frame(dep_flt@test_result)
                sig$Gene_ID <- rownames(dep_flt@test_result)
                sig <- dplyr::left_join(sig, gene_map, by = "Gene_ID")
                sig_genes <- sig$Gene_Name[sig$significant]
                df_filtered <- df[df$Gene_Name %in% sig_genes, , drop = FALSE]
            } else {
                rd <- SummarizedExperiment::rowData(dep_flt)
                sig_genes <- rd$Gene_Name[rd$significant]
                df_filtered <- res$df[res$df$Gene_Name %in% sig_genes, , drop = FALSE]
            }

            DT::datatable(df_filtered,
                extensions = "Buttons",
                options = list(
                    scrollX = TRUE, pageLength = 10, dom = "Bfrtip",
                    buttons = c("copy", "csv", "excel", "pdf", "print")
                )
            )
        })

        # Show optimal k (computed from same filtered matrix)
        output$optimal_k <- renderUI({
            res <- get_dep_result()
            req(res)
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)
            tags$ul(tags$li(paste("Suggested k (from filtered matrix):", res$optimal_k)))
        })
    })
}
