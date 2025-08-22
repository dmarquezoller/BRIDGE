#' @export
processed_integration <- function(input, output, session, rv) {
    shiny::observeEvent(input$process_integrate_data, {
        shiny::req(input$processed_integration)

        selected_tables <- input$processed_integration
        filtered_genes <- list()
        dim_info <- list()

        for (tbl in selected_tables) {
            contrast <- input[[paste0("pi_comparison_selected_", tbl)]]
            dep <- rv$dep_output[[tbl]]

            # Construct column names for contrast
            lfc_col <- paste0(contrast, "_diff")
            padj_col <- paste0(contrast, "_p.adj")

            # Extract the necessary columns from SummarizedExperiment::rowData
            res <- as.data.frame(SummarizedExperiment::rowData(dep))

            if (rv$datatype[[tbl]] == "rnaseq") {
                mapper <- rv$tables[[tbl]]
                mapper <- mapper[, c("Gene_ID", "Gene_Name")]
                res$Gene_ID <- rownames(res)
                res <- merge(res, mapper, by = "Gene_ID", all.x = TRUE)
            }

            if (!all(c("Gene_Name", lfc_col, padj_col) %in% colnames(res))) {
                shiny::showNotification(paste("Missing columns in table:", tbl), type = "error")
                next
            }

            df <- res[, c("Gene_Name", lfc_col, padj_col)]

            colnames(df) <- c("Gene_Name", "diff", "p.adj") # Standardize column names
            df <- df[!is.na(df$Gene_Name), ]
            df <- df[!is.na(df$diff), ]
            df <- df[!is.na(df$p.adj), ]
            # Apply threshold filtering
            keep <- abs(df$diff) >= input$lfc_thresh_pi & df$p.adj <= input$pval_thresh_pi
            filtered_genes[[tbl]] <- df$Gene_Name[keep]

            original_dim <- dim(SummarizedExperiment::assay(dep))
            filtered_dim <- length(filtered_genes[[tbl]])

            dim_info[[tbl]] <- list(
                original = original_dim,
                filtered = c(filtered_dim, original_dim[2])
            )
        }

        # Intersect gene names across all filtered lists
        all_ids <- lapply(filtered_genes, unlist)
        common_ids <- Reduce(intersect, all_ids)
        if (length(common_ids) == 0) {
            shiny::showNotification("No intersected significant genes found across selected datasets.", type = "error")
            rv$intersected_tables_processed <- NULL
            rv$integration_preview_dims <- NULL
            return()
        }

        # Subset SummarizedExperiment::assays by intersected gene names
        intersected_list <- lapply(selected_tables, function(tbl) {
            data <- rv$tables[[tbl]]
            timepoints <- rv$time_cols[[tbl]]

            mat <- as.matrix(data[timepoints])

            # Generate unique IDs
            if ("pepG" %in% colnames(data)) {
                data$unique_id <- paste(data$Gene_Name, data$pepG, sep = "_")
            } else {
                data$unique_id <- data$Gene_Name
            }

            rownames(mat) <- data$unique_id

            # Subset by intersected gene names (gene_name, not unique_id)
            keep <- data$Gene_Name %in% common_ids
            mat[keep, , drop = FALSE]
        })

        names(intersected_list) <- selected_tables

        # Update dimension info
        for (tbl in selected_tables) {
            dim_info[[tbl]]$intersected <- dim(intersected_list[[tbl]])
            data_for_elbow <- safe_row_scale(intersected_list[[tbl]])
        }

        selected_tables <- names(intersected_list)

        scatter_data <- data.frame(Gene_Name = rownames(intersected_list[[1]]))

        for (tbl in selected_tables) {
            contrast <- input[[paste0("pi_comparison_selected_", tbl)]]
            dep <- rv$dep_output[[tbl]]
            lfc_col <- paste0(contrast, "_diff")

            res <- as.data.frame(SummarizedExperiment::rowData(dep))

            if (rv$datatype[[tbl]] == "rnaseq") {
                mapper <- rv$tables[[tbl]]
                mapper <- mapper[, c("Gene_ID", "Gene_Name")]
                res$Gene_ID <- rownames(res)
                res <- merge(res, mapper, by = "Gene_ID", all.x = TRUE)
            }

            df <- res[, c("Gene_Name", lfc_col)]
            colnames(df) <- c("Gene_Name", paste0("LFC: ", tbl))
            scatter_data <- merge(scatter_data, df, by = "Gene_Name", all.x = TRUE)
        }


        lfc_cols <- setdiff(colnames(scatter_data), "Gene_Name")

        plot_list <- list()
        for (i in 1:(length(lfc_cols) - 1)) {
            for (j in (i + 1):length(lfc_cols)) {
                df_plot <- scatter_data %>%
                    select(Gene_Name, x = all_of(lfc_cols[i]), y = all_of(lfc_cols[j])) %>%
                        filter(!is.na(x) & !is.na(y))

                p <- ggplot(df_plot, aes(x = x, y = y, text = Gene_Name)) +
                    geom_point(alpha = 0.7, color = "#2b8cbe") +
                        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
                        labs(
                            title = paste(lfc_cols[i], "vs", lfc_cols[j]),
                            x = lfc_cols[i],
                            y = lfc_cols[j]
                        ) +
                        theme_minimal()

                # Store plotly version
                plot_list[[paste(lfc_cols[i], lfc_cols[j], sep = "_vs_")]] <- p
            }
        }


        # Save results
        rv$scatter_plots <- plot_list
        rv$intersected_tables_processed <- intersected_list
        rv$integration_preview_dims <- lapply(selected_tables, function(tbl) {
            list(
                original = dim(SummarizedExperiment::assay(rv$dep_output[[tbl]])),
                filtered = dim_info[[tbl]]$filtered,
                intersected = dim_info[[tbl]]$intersected
            )
        })
        optimal_k <- safe_nbclust(data_for_elbow, k_min = 2, k_max = 10)
        if (is.na(optimal_k)) {
            showNotification("Not enough clean data to estimate k (after filtering).", type = "warning")
            return(invisible())
        }
        rv$optimal_k <- optimal_k
        names(rv$integration_preview_dims) <- selected_tables
    })
}
