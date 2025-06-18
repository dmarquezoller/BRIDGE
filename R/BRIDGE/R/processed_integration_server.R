#' @export
processed_integration <- function(input, output, session, rv){
    shiny::observeEvent(input$process_integrate_data, {
        req(input$processed_integration)

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
            res <- as.data.frame( SummarizedExperiment::rowData(dep))

            if (rv$datatype[[tbl]] == "rnaseq"){
            mapper <- rv$tables[[tbl]]
            mapper <- mapper[, c("Gene_ID", "Gene_Name")]
            res$Gene_ID <- rownames(res)
            res <- merge(res, mapper, by = "Gene_ID", all.x = TRUE)

            }

            if (!all(c("Gene_Name", lfc_col, padj_col) %in% colnames(res))) {
            showNotification(paste("Missing columns in table:", tbl), type = "error")
            next
            }

            df <- res[, c("Gene_Name", lfc_col, padj_col)]

            colnames(df) <- c("Gene_Name", "diff", "p.adj")  # Standardize column names
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
            showNotification("No intersected significant genes found across selected datasets.", type = "error")
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
            data_for_elbow <- intersected_list[[tbl]]
        }

        # Save results
        rv$intersected_tables_processed <- intersected_list
        rv$integration_preview_dims <- lapply(selected_tables, function(tbl) {
            list(
            original = dim(SummarizedExperiment::assay(rv$dep_output[[tbl]])),
            filtered = dim_info[[tbl]]$filtered,
            intersected = dim_info[[tbl]]$intersected
            )
        })
        elbow <- NbClust(data_for_elbow, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")
        optimal_k <- as.numeric(names(sort(table(elbow$Best.nc[1, ]), decreasing = TRUE)[1]))
        rv$optimal_k <- optimal_k
        names(rv$integration_preview_dims) <- selected_tables
    })
}