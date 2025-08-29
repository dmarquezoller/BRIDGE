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
            sig <- SummarizedExperiment::assay(sig)
            mat <- as.matrix(sig)
            # message("MAT:", head(mat))
            # sanitize matrix (rows/cols need finite variance and ≥2 finite values)
            row_ok <- rowSums(is.finite(mat)) >= 2 & apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
            col_ok <- colSums(is.finite(mat)) >= 2 & apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
            mat <- mat[row_ok, col_ok, drop = FALSE]
            mat_scaled <- safe_row_scale(mat)
            # message("MATFILT:", head(mat_scaled))
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

                    # table snapshot based on same filtered/significant set
                    sig <- DEP2::get_signicant(dep_flt) # only significant rows
                    gene_info <- as.data.frame(SummarizedExperiment::rowData(sig))
                    df <- cbind(gene_info, as.data.frame(mat))
                    df <- df[, c(colnames(gene_info), colnames(mat)), drop = FALSE]
                    df <- stats::na.omit(df)
                    # message("Finishing heatmap task")
                    list(optimal_k = optimal_k, df = df)
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
                heatmap_ready(TRUE)

                ## Now we proceed
                # res <- get_dep_result()
                # req(res)
                # params <- req(last_params())
                ## message("Got DEP result")
                # if (!is.null(key) && !cache$exists(key)) cache$set(key, res)
                #
                # dep_flt_list <- get_depflt(params)
                # dep_flt <- dep_flt_list$dep_flt
                # if (params$clustering) {
                #    k <- as.integer(params$k)
                #    k_max <- max(2L, nrow(dep_flt) - 1L)
                #    if (!is.finite(k) || k < 2L) k <- 2L
                #    if (k > k_max) {
                #        showNotification(sprintf(
                #            "Reducing k from %d to %d (only %d rows).",
                #            input$num_clusters, k_max, nrow(dep_flt)
                #        ), type = "warning")
                #        k <- k_max
                #    }
                # }

                ## consistent y-limits across clusters
                # mat <- dep_flt_list$mat_scaled
                # y_lim <- range(mat, na.rm = TRUE)
                # if (!is.null(k)) {
                #    ht <- ComplexHeatmap::Heatmap(
                #        mat,
                #        row_km = k,
                #        show_row_names = FALSE,
                #        cluster_columns = TRUE
                #    )
                #
                #    # custom per-cluster profile panel (grey lines per gene + thick median)
                #    profiles <- ComplexHeatmap::rowAnnotation(
                #        profile = function(index) {
                #            # index are the row indices for this cluster slice
                #            sub <- mat[index, , drop = FALSE]
                #            n <- ncol(sub)
                #
                #            # set plotting scales: x = columns (samples), y = expression range
                #            pushViewport(viewport(xscale = c(1, n), yscale = y_lim))
                #
                #            # light-grey line per gene
                #            apply(sub, 1, function(v) {
                #                grid.lines(
                #                    x = seq_len(n), y = v, default.unit = "native",
                #                    gp = gpar(col = "#00000022", lwd = 0.6)
                #                ) # translucent grey
                #            })
                #
                #            # overlay median profile (or use colMeans(sub))
                #            med <- apply(sub, 2, median, na.rm = TRUE)
                #            grid.lines(
                #                x = seq_len(n), y = med, default.unit = "native",
                #                gp = gpar(col = "black", lwd = 2)
                #            )
                #
                #            # (optional) axes box
                #            grid.rect()
                #            popViewport()
                #        },
                #        width = unit(40, "mm")
                #    )
                #    ht <- ht + profiles
                # } else {
                #    ht <- ComplexHeatmap::Heatmap(
                #        mat,
                #        cluster_rows = FALSE,
                #        show_row_names = FALSE,
                #        cluster_columns = TRUE
                #    )
                # }
                # ht
                # ht <- ComplexHeatmap::draw(ht, merge_legend = TRUE, newpage = FALSE)
                # roword(ComplexHeatmap::row_order(ht))
                # ht
                # heatmap_ready(TRUE)

                # session$onFlushed(function() {
                #    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
                #        input, output, session,
                #        ht_list = ht,
                #        heatmap_id = ich_id,
                #        click_action = click_action,
                #        brush_action = brush_action,
                #        hover_action = NULL,
                #        res = 96,
                #        show_cell_fun = TRUE,
                #        show_layer_fun = TRUE
                #    )
                # }, once = TRUE)
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

        output$ht <- renderPlot(
            {
                req(heatmap_ready())
                key <- rv$current_dep_heatmap_key[[tbl_name]]
                # Now we proceed
                res <- get_dep_result()
                req(res)
                params <- req(last_params())
                # message("Got DEP result")
                if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

                dep_flt_list <- get_depflt(params)
                dep_flt <- dep_flt_list$dep_flt
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
                # consistent y-limits across clusters
                mat <- dep_flt_list$mat_scaled
                local_mat <- mat
                y_lim <- range(local_mat, na.rm = TRUE)
                if (!all(is.finite(y_lim)) || diff(y_lim) == 0) y_lim <- c(-1, 1)

                if (!is.null(k)) {
                    ht <- ComplexHeatmap::Heatmap(
                        local_mat,
                        row_km = k,
                        show_row_names = FALSE,
                        cluster_columns = TRUE
                    ) + ComplexHeatmap::rowAnnotation(
                        profile = ComplexHeatmap::anno_empty(which = "row", width = grid::unit(50, "mm")),
                        annotation_name_gp = grid::gpar(col = NA)
                    )
                } else {
                    ht <- ComplexHeatmap::Heatmap(
                        local_mat,
                        cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = TRUE
                    )
                }

                # custom per-cluster profile panel (grey lines per gene + thick median)
                # profiles <- ComplexHeatmap::rowAnnotation(
                #    profile = ComplexHeatmap::anno_lines(
                #        mat,  # rows=genes, cols=samples
                #        ylim = y_lim,
                #        gp   = grid::gpar(lwd = 0.6),
                #        pch  = NA,                       # no points
                #        axis = TRUE
                #    ),
                #    width = grid::unit(40, "mm")
                # )

                ht <- ComplexHeatmap::draw(ht, merge_legend = TRUE, newpage = TRUE)
                roword(ComplexHeatmap::row_order(ht))

                y_lim <- range(local_mat, na.rm = TRUE)
                if (!all(is.finite(y_lim)) || diff(y_lim) == 0) y_lim <- c(-1, 1)

                rord <- ComplexHeatmap::row_order(ht)
                cord <- ComplexHeatmap::column_order(ht)
                labs <- colnames(local_mat)[cord]
                xs <- seq_along(cord)

                for (s in seq_along(rord)) {
                    idx <- rord[[s]]
                    sub <- local_mat[idx, cord, drop = FALSE]
                    if (!is.matrix(sub) || nrow(sub) == 0L) next

                    med <- apply(sub, 2, stats::median, na.rm = TRUE)
                    q1 <- apply(sub, 2, stats::quantile, probs = 0.25, na.rm = TRUE)
                    q3 <- apply(sub, 2, stats::quantile, probs = 0.75, na.rm = TRUE)

                    # message("ANNO: ", idx, head(med), head(q1), head(q3))

                    ComplexHeatmap::decorate_annotation("profile", slice = s, {
                        # 1) enter the scaled, clipped viewport
                        grid::pushViewport(grid::viewport(
                            xscale = c(1, length(xs)),
                            yscale = y_lim,
                            clip   = "on"
                        ))
                        on.exit(grid::popViewport(), add = TRUE)

                        # (optional) red border in the scaled space (so you see clipping is active)
                        grid::grid.rect(gp = grid::gpar(fill = NA, col = "red"))

                        # 2) sanity line IN 'native' coords — should be mid (y=0) if y_lim spans negatives
                        grid::grid.lines(
                            x = grid::unit(c(1, length(xs)), "native"),
                            y = grid::unit(c(0, 0), "native"),
                            gp = grid::gpar(col = "blue")
                        )

                        # 3) ribbon (only finite columns)
                        ok <- is.finite(q1) & is.finite(q3)
                        if (any(ok)) {
                            x_ok <- xs[ok]
                            grid::grid.polygon(
                                x = grid::unit(c(x_ok, rev(x_ok)), "native"),
                                y = grid::unit(c(q1[ok], rev(q3[ok])), "native"),
                                gp = grid::gpar(fill = "#00000022", col = NA)
                            )
                        }

                        # 4) median trend (draw in finite runs)
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

                        # 5) ticks + rotated labels (still inside same viewport)
                        grid::grid.xaxis(at = xs, label = FALSE)
                        for (i in xs) {
                            grid::grid.text(
                                labs[i],
                                x = grid::unit(i, "native"),
                                y = grid::unit(0, "npc") + grid::unit(1.2, "mm"),
                                rot = 90, just = c(1, 0.5), gp = grid::gpar(cex = 0.55)
                            )
                        }
                    })
                }
                invisible(NULL)
            },
            res = 96
        )


        output$ht_sig <- DT::renderDT({
            params <- req(last_params())
            roword <- req(roword())
            res <- get_dep_result()
            req(res)
            key <- rv$current_dep_heatmap_key[[tbl_name]]
            if (!is.null(key) && !cache$exists(key)) cache$set(key, res)
            cluster.all <- NULL
            if (!is.null(roword)) {
                cluster.all <- list()
                # loop to extract genes for each cluster.
                for (i in 1:length(roword)) {
                    # message("clusters: ", i, head(roword[[i]]))
                    if (i == 1) {
                        clu <- t(t(res$df[roword[[i]], ]$Gene_ID))
                        cluster.all <- cbind(clu, paste("cluster", i, sep = ""))
                    } else {
                        clu <- t(t(res$df[roword[[i]], ]$Gene_ID))
                        clu <- cbind(clu, paste("cluster", i, sep = ""))
                        cluster.all <- rbind(cluster.all, clu)
                    }
                }
                colnames(cluster.all) <- c("Gene_ID", "Cluster")
                cluster.all <- as.data.frame(cluster.all)
            }

            # message("clusters: ", head(cluster.all))
            dep_flt_list <- get_depflt(params)
            dep_flt <- dep_flt_list$dep_flt
            if (methods::is(dep_flt, "DEGdata")) {
                df <- res$df
                df$name <- rownames(df)
                df$Gene_ID <- gsub("^(.*)_", "", rownames(df))
                gene_map <- rv$tables[[tbl_name]][, c("Gene_ID", "Gene_Name")]
                df <- dplyr::left_join(df, gene_map, by = "Gene_ID")
                if (!is.null(cluster.all)) {
                    df <- dplyr::left_join(df, cluster.all, by = "Gene_ID")
                }

                sig <- as.data.frame(dep_flt@test_result)
                sig$Gene_ID <- gsub("^(.*)_", "", rownames(dep_flt@test_result))
                sig <- dplyr::left_join(sig, gene_map, by = "Gene_ID")
                sig_genes <- sig$Gene_Name[sig$significant]

                df_filtered <- df[df$Gene_Name %in% sig_genes, , drop = FALSE]
                rownames(df_filtered) <- NULL
                if (!is.null(cluster.all)) {
                    df_filtered <- df_filtered %>%
                        select(Gene_ID, Gene_Name, Cluster, everything()) %>%
                        mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                } else {
                    df_filtered <- df_filtered %>%
                        select(Gene_ID, Gene_Name, everything()) %>%
                        mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                }
            } else {
                rd <- SummarizedExperiment::rowData(dep_flt)
                sig_genes <- rd$Gene_Name[rd$significant]
                df_filtered <- res$df[res$df$Gene_Name %in% sig_genes, , drop = FALSE] %>% dplyr::select(-ID)
                if (!is.null(cluster.all)) {
                    df_filtered <- dplyr::left_join(df_filtered, cluster.all, by = "Gene_ID")
                }
                rownames(df_filtered) <- NULL
                if (!is.null(cluster.all)) {
                    df_filtered <- df_filtered %>%
                        select(Gene_ID, Gene_Name, Cluster, everything()) %>%
                        mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                } else {
                    df_filtered <- df_filtered %>%
                        select(Gene_ID, Gene_Name, everything()) %>%
                        mutate(Gene_Name = stringr::str_to_title(Gene_Name))
                }
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
