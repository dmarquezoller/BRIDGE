#' @export
int_heatmap_server <- function(input, output, session, rv) {
    output$integrated_heatmaps <- shiny::renderUI({
        shiny::req(rv$intersected_matrix_processed, input$heatmap_k)
        tables <- rv$intersected_matrix_processed
        htmltools::tagList(
            lapply(names(rv$intersected_matrix_processed), function(tbl) {
                plotOutput(outputId = paste0("heatmap_", tbl), width = "80%")
            }),
            shiny::br(),
            lapply(names(tables), function(tbl) {
                shinydashboard::box(
                    title = paste("Clusters Table:", tbl),
                    width = 12,
                    solidHeader = TRUE,
                    status = "info",
                    collapsible = TRUE,
                    collapsed = FALSE,
                    DT::DTOutput(outputId = paste0("cluster_table_", tbl))
                )
            })
        )
    })

    shiny::observe({
        shiny::req(rv$intersected_matrix_processed, input$heatmap_k)

        all_tables <- rv$intersected_matrix_processed

        first_tbl <- all_tables[[1]]
        mat_scaled <- safe_row_scale(first_tbl) # Scale across rows
        # allow some NAs; require ≥2 finite values and non-zero SD
        row_ok <- rowSums(is.finite(mat_scaled)) >= 2 &
            apply(mat_scaled, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)

        mat_scaled <- mat_scaled[row_ok, , drop = FALSE]
        if (nrow(mat_scaled) < 2 || ncol(mat_scaled) < 2) {
            shiny::showNotification("Too few usable rows/columns after cleaning.", type = "warning")
            return(NULL)
        }

        k <- as.integer(input$heatmap_k)
        if (!is.finite(k) || k < 2) k <- 2

        # maximum valid k is nrow(mat) - 1
        k_max <- max(2L, nrow(mat_scaled) - 1L)
        if (k > k_max) {
            shiny::showNotification(
                sprintf("Reducing k from %d to %d (only %d rows after filtering).", k, k_max, nrow(mat_scaled)),
                type = "warning"
            )
            k <- k_max
        }

        # still not enough rows?
        if (k_max < 2L) {
            shiny::showNotification("Not enough rows for clustering (need ≥ 2).", type = "error")
            return(NULL)
        }

        # safe kmeans
        km <- try(stats::kmeans(mat_scaled, centers = k, nstart = 10, iter.max = 100), silent = TRUE)
        if (inherits(km, "try-error")) {
            shiny::showNotification("k-means failed on cleaned matrix; try smaller k or relax filters.", type = "error")
            return(NULL)
        }

        # Get cluster labels and order of genes
        cluster_labels <- km$cluster
        ordered_genes <- rownames(mat_scaled)[order(cluster_labels)]

        # Save info for use across plots
        rv$heatmap_clusters <- list(
            order = ordered_genes,
            cluster = cluster_labels
        )

        gene_map <- data.frame(unique_id = rownames(first_tbl), Gene_Name = sapply(strsplit(rownames(first_tbl), "_"), `[`, 1))
        rv$gene_to_cluster <- setNames(rv$heatmap_clusters$cluster, gene_map$Gene_Name)

        # Render heatmaps
        lapply(names(all_tables), function(tbl) {
            local({
                tbl_name <- tbl
                mat <- all_tables[[tbl_name]]
                mat_scaled_tbl <- safe_row_scale(mat)

                gene_names <- sapply(strsplit(rownames(mat_scaled_tbl), "_"), `[`, 1)
                cluster_vec <- factor(rv$gene_to_cluster[gene_names], levels = 1:k)

                valid <- !is.na(cluster_vec)
                mat_ordered <- mat_scaled_tbl[valid, , drop = FALSE]
                cluster_vec <- cluster_vec[valid]

                ordered <- order(cluster_vec)
                mat_ordered <- mat_ordered[ordered, , drop = FALSE]
                cluster_vec <- cluster_vec[ordered]

                output[[paste0("heatmap_", tbl_name)]] <- shiny::renderPlot({
                    mat <- all_tables[[tbl_name]]
                    mat_scaled_tbl <- safe_row_scale(mat)

                    # Build cluster average profiles
                    cluster_ids <- split(seq_len(nrow(mat_ordered)), cluster_vec)
                    line_profiles <- t(vapply(cluster_ids, function(idxs) {
                        colMeans(mat_ordered[idxs, , drop = FALSE], na.rm = TRUE)
                    }, FUN.VALUE = numeric(ncol(mat_ordered))))

                    # Safeguard: avoid empty trend lines
                    if (nrow(line_profiles) == 0 || ncol(line_profiles) == 0) {
                        shiny::showNotification("Could not compute trend lines — data missing or clustering failed", type = "error")
                        return(NULL)
                    }

                    # Normalize trend lines for plotting
                    line_profiles_norm <- t(apply(line_profiles, 1, function(x) {
                        rng <- range(x, na.rm = TRUE)
                        if (diff(rng) == 0) rep(0.5, length(x)) else (x - rng[1]) / diff(rng)
                    }))

                    # Annotation: line plots next to each cluster
                    trend_anno <- ComplexHeatmap::rowAnnotation(trend = ComplexHeatmap::anno_link(
                        align_to = cluster_vec,
                        which = "row",
                        panel_fun = function(index, nm) {
                            grid::grid.rect()
                            grid::grid.lines(
                                x = seq_len(ncol(line_profiles_norm)) / ncol(line_profiles_norm),
                                y = line_profiles_norm[as.integer(nm), ],
                                gp = ggfun::gpar(col = "#2b8cbe", lwd = 1)
                            )
                        },
                        side = "right",
                        size = unit(3, "cm"),
                        width = unit(5, "cm")
                    ))
                    # Final heatmap with trend annotation
                    ComplexHeatmap::Heatmap(
                        mat_ordered,
                        name = tbl_name,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_dend = FALSE,
                        show_column_dend = TRUE,
                        show_row_names = FALSE,
                        row_split = cluster_vec,
                        row_names_gp = ggfun::gpar(fontsize = 6),
                        column_names_gp = ggfun::gpar(fontsize = 8),
                        heatmap_legend_param = list(title = "Expression"),
                        right_annotation = trend_anno
                    )
                })
            })
        })

        lapply(names(all_tables), function(tbl) {
            local({
                tbl_name <- tbl
                data <- all_tables[[tbl_name]]

                gene_names <- sapply(strsplit(rownames(data), "_"), `[`, 1)
                clusters <- rv$gene_to_cluster[gene_names]
                clusters[is.na(clusters)] <- NA

                df <- data.frame(
                    Unique_ID = rownames(data),
                    Gene_Name = gene_names,
                    Cluster = clusters,
                    stringsAsFactors = FALSE
                )

                output[[paste0("cluster_table_", tbl_name)]] <- DT::renderDT({
                    DT::datatable(df, extensions = "Buttons", options = list(scrollX = TRUE, pageLength = 5, dom = "Bfrtip", buttons = c("copy", "csv", "excel", "pdf", "print")))
                })
            })
        })


        output$lfc_scatter_selector <- shiny::renderUI({
            shiny::req(rv$scatter_plots)
            shiny::selectInput("scatter_comparisons", "Select comparison:", choices = names(rv$scatter_plots))
        })

        output$lfc_scatter_ui <- shiny::renderUI({
            shiny::req(rv$intersected_tables_processed)
            selected_tables <- names(rv$intersected_tables_processed)
            shiny::req(rv$scatter_plots)

            if (any(rv$datatype[selected_tables] == "phosphoproteomics")) {
                # Return a nicely styled message div
                shiny::div(
                    style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                    "Scatter plot is not available when phosphoproteomics data is included because of non-1:1 row mapping."
                )
            } else {
                plotly::plotlyOutput("lfc_scatter_plot", height = "400px", width = "100%")
            }
        })



        output$lfc_scatter_plot <- plotly::renderPlotly({
            shiny::req(rv$scatter_plots)
            shiny::req(input[["scatter_comparisons"]])
            plotly::ggplotly(rv$scatter_plots[[input[["scatter_comparisons"]]]], tooltip = "text")
        })
    })
}
