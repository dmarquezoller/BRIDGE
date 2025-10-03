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
                datapoint_cols <- raw_input$datapoint_cols
                datatype <- raw_input$datatype

                clean <- as.matrix(raw_data[datapoint_cols])
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
                    datapoint_cols = isolate(rv$data_cols[[tbl_name]]),
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
        last_params <- reactiveVal(NULL) # frozen on click
        last_mat <- reactiveVal(NULL)
        roword <- reactiveVal(NULL)
        ns <- session$ns
        ich_id <- "dep_ht"

        # InteractiveComplexHeatmapAction click
        click_action <- function(df, output) {
            # df has heatmap name, row_index, column_index, etc.
            output[["cell_info"]] <- renderUI({
                if (is.null(df)) {
                    return(NULL)
                }
                HTML(sprintf(
                    "<b>Clicked:</b> heatmap <code>%s</code> — row %s, column %s",
                    df$heatmap, df$row_index, df$column_index
                ))
            })
        }
        # InteractiveComplexHeatmapAction brush
        brush_action <- function(df, output) {
            r <- unique(unlist(df$row_index))
            c <- unique(unlist(df$column_index))
            output[["selection_table"]] <- renderUI({
                m <- last_mat()
                if (is.null(m) || length(r) == 0 || length(c) == 0) {
                    return(NULL)
                }
                sel <- round(m[r, c, drop = FALSE], 3)
                htmltools::HTML(paste0("<pre>", paste(capture.output(print(sel)), collapse = "\n"), "</pre>"))
            })
        }

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
            # message("Getting DEP result")
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && cache$exists(key)) {
                return(cache$get(key))
            }
            heatmap_task$result()
        }
        get_depflt <- function(params) {
            # message("Getting DEP filtered result")
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
            dep_output <- strip_sig(dep_output)
            dep_flt <- DEP2::add_rejections(dep_output, alpha = params$p_cut, lfc = params$lfc_cut)

            mats <- get_depflt_matrix(dep_flt)
            # message("MATS:", head(mats))
            mat <- mats$mat
            mat_scaled <- mats$mat_scaled
            last_mat(mat_scaled)
            ret_list <- list(mat = mat, mat_scaled = mat_scaled, dep_flt = dep_flt)
            depflt_cache[[key]] <- ret_list
            # message("DEP filtered result cached", str(ret_list))
            ret_list
        }

        get_depflt_matrix <- function(dep_obj) {
            # message("Getting DEP filtered matrix")
            sig <- DEP2::get_signicant(dep_obj) # only significant rows
            # If there are no significant rows, or the object is empty, exit early
            if (is.null(sig) || nrow(sig) == 0L) {
                last_mat(NULL)
                return(list(mat = NULL, mat_scaled = NULL))
            }
            sig <- SummarizedExperiment::assay(sig)
            # If assay is NULL or empty, exit early
            mat <- as.matrix(sig)
            if (is.null(mat) || length(mat) == 0L) {
                last_mat(NULL)
                return(list(mat = NULL, mat_scaled = NULL))
            }
            # message("MAT:", head(mat, n = 3))
            # sanitize matrix (rows/cols need finite variance and ≥2 finite values)
            row_ok <- rowSums(is.finite(mat)) >= 2 & apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
            col_ok <- colSums(is.finite(mat)) >= 2 & apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
            mat <- mat[row_ok, col_ok, drop = FALSE]
            mat_scaled <- safe_row_scale(mat)
            # message("MATFILT:", head(mat_scaled, n = 3))
            last_mat(mat_scaled)
            list(mat = mat, mat_scaled = mat_scaled)
        }

        # Background task: compute optimal_k from the filtered (significant) matrix + a DF snapshot
        heatmap_task <- ExtendedTask$new(function(args) {
            promises::future_promise(
                {
                    dep_output <- args$dep_output
                    p_cut <- args$p_cut
                    lfc_cut <- args$lfc_cut
                    datatype <- args$datatype

                    if (methods::is(dep_output, "DEGdata")) {
                        # RNA-seq case: clean significant columns from test_result
                        tr <- dep_output@test_result
                        tr <- tr[, !grepl("_significant$|^significant$", colnames(tr))]
                        dep_output@test_result <- tr
                    } else {
                        # Proteomics case: clean rowData
                        rd <- SummarizedExperiment::rowData(dep_output)
                        rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
                        SummarizedExperiment::rowData(dep_output) <- rd
                    }

                    # filter with current cutoffs
                    dep_flt <- DEP2::add_rejections(dep_output, alpha = p_cut, lfc = lfc_cut)
                    mats <- get_depflt_matrix(dep_flt)
                    mat <- mats$mat
                    mat_scaled <- mats$mat_scaled
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

                    df <- cbind(as.data.frame(SummarizedExperiment::assay(dep_flt)), as.data.frame(SummarizedExperiment::rowData(dep_flt)))
                    df_table <- df

                    names <- switch(datatype,
                        proteomics        = paste0(stringr::str_to_title(df$Gene_Name), "_", df$Protein_ID),
                        phosphoproteomics = paste0(stringr::str_to_title(df$Gene_Name), "_", df$Protein_ID, "_", df$pepG),
                        rnaseq            = rownames(df)
                    )
                    df$name <- names
                    if (datatype == "rnaseq") {
                        df$Gene_Name <- gsub("_.*", "", df$name, perl = TRUE)
                        df$Gene_ID <- gsub(".*_", "", df$name, perl = TRUE)
                    }
                    df <- df %>% dplyr::select(Gene_ID, Gene_Name, name, dplyr::everything())
                    #message("BEF: ", nrow(df), "\nhead: ", head(df))
                    df <- stats::na.omit(df)
                    #message("AFT: ", nrow(df))
                    list(optimal_k = optimal_k, df = df, table = df_table, datatype = datatype)
                },
                seed = TRUE
            )
        })

        draw_heatmap_task <- ExtendedTask$new(function(args) {
            promises::future_promise(
                {
                    k <- args$optimal_k
                    mat <- args$mat

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

                    # consistent y-limits across clusters
                    y_lim <- range(mat, na.rm = TRUE)
                    if (!is.null(k)) {
                        ht <- ComplexHeatmap::Heatmap(
                            mat,
                            row_km = k,
                            row_km_repeats = 1,
                            show_row_names = FALSE,
                            cluster_columns = TRUE
                        )

                        # custom per-cluster profile panel (grey lines per gene + thick median)
                        profiles <- ComplexHeatmap::rowAnnotation(
                            profile = function(index) {
                                # index are the row indices for this cluster slice
                                sub <- mat[index, , drop = FALSE]
                                n <- ncol(sub)

                                # set plotting scales: x = columns (samples), y = expression range
                                pushViewport(viewport(xscale = c(1, n), yscale = y_lim))

                                # light-grey line per gene
                                apply(sub, 1, function(v) {
                                    grid.lines(
                                        x = seq_len(n), y = v, default.unit = "native",
                                        gp = gpar(col = "#00000022", lwd = 0.6)
                                    ) # translucent grey
                                })

                                # overlay median profile (or use colMeans(sub))
                                med <- apply(sub, 2, median, na.rm = TRUE)
                                grid.lines(
                                    x = seq_len(n), y = med, default.unit = "native",
                                    gp = gpar(col = "black", lwd = 2)
                                )

                                # (optional) axes box
                                grid.rect()
                                popViewport()
                            },
                            width = unit(40, "mm")
                        )
                        ht <- ht + profiles
                    } else {
                        ht <- ComplexHeatmap::Heatmap(
                            mat,
                            cluster_rows = FALSE,
                            show_row_names = FALSE,
                            cluster_columns = TRUE
                        )
                    }
                    # message("Finishing heatmap DRAW task", class(ht))
                    return(ht)
                },
                seed = TRUE
            )
        })

        # Click: freeze params and kick off background compute (key includes cutoffs)
        observeEvent(input$recompute_heatmap,
            {
                req(rv$dep_output[[tbl_name]])
                params <- list(
                    p_cut = input$heatmap_pcutoff, # raw FDR 0..1
                    lfc_cut = input$heatmap_fccutoff,
                    clustering = isTRUE(input$clustering),
                    k = input$num_clusters,
                    dep_output = isolate(rv$dep_output[[tbl_name]]),
                    datatype = isolate(rv$datatype[[tbl_name]]),
                    columns_key = paste(isolate(rv$data_cols[[tbl_name]]), collapse = "_"),
                    dummy = NULL
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
                        p_cut      = params$p_cut,
                        lfc_cut    = params$lfc_cut,
                        datatype   = params$datatype
                    ))
                }
                # Pre-warm depflt_cache for this (table, columns, p, lfc)

                depflt_key <- paste(
                    tbl_name, params$columns_key,
                    sprintf("pcut=%.4f", params$p_cut),
                    sprintf("lfc=%.3f", params$lfc_cut),
                    "depflt",
                    sep = "_"
                )

                if (is.null(depflt_cache[[depflt_key]])) {
                    pr <- promises::future_promise(
                        {
                            # Do NOT touch reactiveValues/last_mat in the future thread
                            dep_output <- strip_sig(params$dep_output)
                            dep_flt <- DEP2::add_rejections(dep_output, alpha = params$p_cut, lfc = params$lfc_cut)

                            # Build matrices without using get_depflt_matrix (it updates last_mat)
                            sig <- DEP2::get_signicant(dep_flt)
                            if (is.null(sig) || nrow(sig) == 0L) {
                                return(list(mat = NULL, mat_scaled = NULL, dep_flt = dep_flt))
                            }

                            mat <- as.matrix(SummarizedExperiment::assay(sig))
                            row_ok <- rowSums(is.finite(mat)) >= 2 &
                                apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
                            col_ok <- colSums(is.finite(mat)) >= 2 &
                                apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
                            if (!any(row_ok) || !any(col_ok)) {
                                return(list(mat = NULL, mat_scaled = NULL, dep_flt = dep_flt))
                            }

                            mat <- mat[row_ok, col_ok, drop = FALSE]
                            mat_scaled <- safe_row_scale(mat)

                            list(mat = mat, mat_scaled = mat_scaled, dep_flt = dep_flt)
                        },
                        seed = TRUE
                    )

                    # Commit the result to reactiveValues on the main Shiny thread
                    pr <- promises::then(pr, function(ret) {
                        depflt_cache[[depflt_key]] <- ret
                        last_mat(ret$mat_scaled) # optional: keep preview matrix in sync
                        NULL
                    })

                    # Swallow/log errors so UI doesn’t break
                    pr <- promises::catch(pr, function(e) {
                        # message("Prewarm depflt failed: ", conditionMessage(e))
                        NULL
                    })
                }
                heatmap_ready(TRUE)
            },
            ignoreInit = TRUE
        )

        # Spinner container
        output$ht_slot <- renderUI({
            if (!heatmap_ready()) {
                return(div(
                    style = "padding:10px; color:#777;",
                    "Choose settings and click “Compute heatmap”."
                ))
            }
            shinycssloaders::withSpinner(
                plotOutput(session$ns("ht"), height = "800px", width = "100%"),
                type = 8, color = "#2b8cbe", caption = "Loading..."
            )
            # InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
            #    heatmap_id = ns(ich_id),
            #    width1 = 700, height1 = 520,
            #    width2 = 520, height2 = 380,
            #    layout = "1|(2-3)"
            # )
        })


        output$ht <- renderPlot({
            req(heatmap_ready())
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            res <- get_dep_result()
            req(res)
            params <- req(last_params())
            # message("Got DEP result")
            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

            dep_flt_list <- get_depflt(params)
            dep_flt <- dep_flt_list$dep_flt
# make sure k is defined even when clustering is off
k <- NULL
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
            }
local_mat <- dep_flt_list$mat_scaled
# handle empty/NULL matrix
if (is.null(local_mat) || nrow(local_mat) == 0L || ncol(local_mat) == 0L) {
    plot.new()
    text(0.5, 0.5,
        "No significant features with the current FDR/log2FC cutoffs.",
        cex = 1.1
    )
    return(invisible())
}

            # consistent y-limits across clusters
            y_lim <- range(local_mat, na.rm = TRUE)
            if (!all(is.finite(y_lim)) || diff(y_lim) == 0) y_lim <- c(-1, 1)

            if (!is.null(k)) {
                ht <- ComplexHeatmap::Heatmap(
                    local_mat,
                    row_km = k,
                    show_row_names = FALSE,
                    cluster_columns = FALSE
                ) + ComplexHeatmap::rowAnnotation(
                    profile = ComplexHeatmap::anno_empty(which = "row", width = grid::unit(55, "mm"), height = grid::unit(30, "mm")),
                    annotation_name_gp = grid::gpar(col = NA)
                )
            } else {
                ht <- ComplexHeatmap::Heatmap(
                    local_mat,
                    cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE
                )
            }

            # custom per-cluster profile panel
            ht <- ComplexHeatmap::draw(ht, merge_legend = TRUE, newpage = TRUE)
            rord <- ComplexHeatmap::row_order(ht)
            roword(rord)
            if (is.list(rord)) {
                n_slices <- length(rord)
            } else {
                rord <- list(rord)
                n_slices <- 1L
            }

            cord <- ComplexHeatmap::column_order(ht)
            if (is.list(cord)) cord <- cord[[1L]]

            labs <- colnames(local_mat)[cord]
            xs <- seq_along(cord)

            y_lim <- range(local_mat, na.rm = TRUE)
            if (!all(is.finite(y_lim)) || diff(y_lim) == 0) y_lim <- c(-1, 1)

            for (s in seq_len(n_slices)) {
                idx <- rord[[s]]
                sub <- local_mat[idx, cord, drop = FALSE]
                if (!is.matrix(sub) || nrow(sub) == 0L) next

                med <- apply(sub, 2, stats::median, na.rm = TRUE)
                q1 <- apply(sub, 2, stats::quantile, probs = 0.25, na.rm = TRUE)
                q3 <- apply(sub, 2, stats::quantile, probs = 0.75, na.rm = TRUE)

                ComplexHeatmap::decorate_annotation("profile", slice = s, {
                    # Scaled, clipped viewport for drawing
                    grid::pushViewport(grid::viewport(
                        xscale = c(0.5, length(xs) + 0.5),
                        yscale = y_lim,
                        clip   = "on"
                    ))
                    on.exit(grid::popViewport(), add = TRUE)

                    # Midline at y = 0 (if in range)
                    grid::grid.lines(
                        x = grid::unit(c(1, length(xs)), "native"),
                        y = grid::unit(c(0, 0), "native"),
                        gp = grid::gpar(col = "blue")
                    )

                    # Ribbon (IQR)
                    ok <- is.finite(q1) & is.finite(q3)
                    if (any(ok)) {
                        x_ok <- xs[ok]
                        grid::grid.polygon(
                            x = grid::unit(c(x_ok, rev(x_ok)), "native"),
                            y = grid::unit(c(q1[ok], rev(q3[ok])), "native"),
                            gp = grid::gpar(fill = "#00000022", col = NA)
                        )
                    }

                    # Median trend (draw in finite runs)
                    okm <- is.finite(med)
                    if (any(okm)) {
                        runs <- split(xs[okm], cumsum(c(1, diff(xs[okm]) != 1)))
                        for (r in runs) {
                            grid::grid.lines(
                                x = grid::unit(r, "native"),
                                y = grid::unit(med[r], "native"),
                                gp = grid::gpar(col = "black", lwd = 2)
                            )
                        }
                        grid::grid.points(
                            x = grid::unit(xs[okm], "native"),
                            y = grid::unit(med[okm], "native"),
                            pch = 16, size = grid::unit(0.7, "mm")
                        )
                    } else {
                        grid::grid.text("no data", gp = grid::gpar(cex = 0.6, col = "grey50"))
                    }

                    # Only draw ticks + labels for the bottom-most slice
                    if (s == n_slices) {
                        grid::grid.xaxis(at = xs, label = FALSE)

                        for (i in xs) {
                            # measure THIS label at rot=90, cex=0.55
                            tg <- grid::textGrob(labs[i], gp = grid::gpar(cex = 0.55), rot = 90)
                            w_nat <- grid::convertX(grid::grobWidth(tg), "native") # used to nudge edges
                            h_mm <- grid::convertY(grid::grobHeight(tg), "mm", valueOnly = TRUE) # for vertical clearance

                            # default x at the tick
                            x_pos <- grid::unit(i, "native")
                            if (i == xs[1]) {
                                x_pos <- x_pos + 0.2 * w_nat # nudge first label inside
                            } else if (i == xs[length(xs)]) {
                                x_pos <- x_pos - 0.2 * w_nat # nudge last label inside
                            }

                            # center of the label sits at half its height + margin above the bottom
                            y_pos <- grid::unit(h_mm / 2 + 1.2, "mm") # increase 1.2 -> 2 mm if needed

                            grid::grid.text(
                                labs[i],
                                x = x_pos,
                                y = y_pos,
                                rot = 90,
                                just = c(0.5, 0.5), # position by the label's center
                                gp = grid::gpar(cex = 0.55)
                            )
                        }
                    }
                })
                #
            }
        })


        output$ht_panel <- renderUI({
            DT::DTOutput(session$ns("ht_sig"), height = "300px")
        })

        output$ht_sig <- DT::renderDT({
            # message("Fetching DEP result for sig table")
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            res <- req(get_dep_result())
            # message("Got DEP result for sig table")
            params <- req(last_params())
            roword <- req(roword())
            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)
            # message("Getting DEP flt result for sig table")
            dep_flt_list <- get_depflt(params)
            dep_flt <- dep_flt_list$dep_flt
            # message("Got DEP flt result for sig table")
            cluster.all <- NULL
            if (!is.null(roword)) {
                cluster.all <- list()
                # loop to extract genes for each cluster.
                for (i in 1:length(roword)) {
                    if (i == 1) {
                        clu <- t(t(res$df[roword[[i]], ]$name))
                        cluster.all <- cbind(clu, paste("cluster", i, sep = ""))
                    } else {
                        clu <- t(t(res$df[roword[[i]], ]$name))
                        clu <- cbind(clu, paste("cluster", i, sep = ""))
                        cluster.all <- rbind(cluster.all, clu)
                    }
                }
                colnames(cluster.all) <- c("name", "Cluster")
                cluster.all <- as.data.frame(cluster.all)
            }
            
            df <- res$df
            # message("DF: ", str(res$df))        
            if (methods::is(dep_flt, "DEGdata")) {
                if (!is.null(cluster.all)) {
                    df <- dplyr::left_join(df, cluster.all, by = "name")
                }
                sig <- as.data.frame(dep_flt@test_result)
                sig$Gene_ID <- gsub("^.*_", "", rownames(sig))
                sig$Gene_Name <- gsub("_.*$", "", rownames(sig))
                sig_genes <- sig$Gene_Name[sig$significant]
                df_filtered <- df[stringr::str_to_lower(df$Gene_Name) %in% stringr::str_to_lower(sig_genes), , drop = FALSE]
                rownames(df_filtered) <- NULL
                if (!is.null(cluster.all)) {
                    df_filtered <- df_filtered %>%
                        dplyr::select(Gene_ID, Gene_Name, name, Cluster, everything()) %>%
                            dplyr::mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                } else {
                    df_filtered <- df_filtered %>%
                        dplyr::select(Gene_ID, Gene_Name, name, everything()) %>%
                            dplyr::mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                }
            } else {
                rd <- SummarizedExperiment::rowData(dep_flt)
            
                sig_genes <- rd$Gene_Name[rd$significant]
                df_filtered <- res$df[stringr::str_to_lower(res$df$Gene_Name) %in% stringr::str_to_lower(sig_genes), , drop = FALSE] %>% dplyr::select(-ID)
                if (!is.null(cluster.all)) {
                    df_filtered <- dplyr::left_join(df_filtered, cluster.all, by = "name")
                }
                rownames(df_filtered) <- NULL
                if (!is.null(cluster.all)) {
                    df_filtered <- df_filtered %>%
                        dplyr::select(Gene_ID, Gene_Name, name, Cluster, everything()) %>%
                        dplyr::mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                } else {
                    df_filtered <- df_filtered %>%
                        dplyr::select(Gene_ID, Gene_Name, name, everything()) %>%
                        dplyr::mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                }
            }
            # message("DF filtered: ", nrow(df_filtered), " rows")
            DT::datatable(df_filtered %>% dplyr::select(where(~!is.numeric(.)), where(is.numeric)),
                extensions = "Buttons",
                filter = "top",
                options = list(scrollX = TRUE, processing = TRUE, pageLength = 10, dom = "Bfrtip", buttons = c("copy", "csv", "excel", "pdf", "print"))
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
