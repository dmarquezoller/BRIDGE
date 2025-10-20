#' @export
processed_integration <- function(input, output, session, rv) {
    observeEvent(input$process_integrate_data, {
        req(input$processed_integration)

        selected_tables <- input$processed_integration
        filtered_ids <- list()
        dim_info <- list()

        extract_features <- function(dep, datatype, contrast) {
            rd <- as.data.frame(SummarizedExperiment::rowData(dep))
            # Add Gene_Name/Gene_ID consistently
            if (datatype == "rnaseq") {
                rd <- as.data.frame(dep@test_result)
                rd$Gene_ID <- gsub("^.*_", "", rownames(rd))
                rd$Gene_Name <- gsub("_.*$", "", rownames(rd))
            } else {
                rd <- as.data.frame(SummarizedExperiment::rowData(dep))
            }
            rownames(rd) <- NULL
            rd <- rd %>%
                dplyr::select(where(~ !is.numeric(.)), where(is.numeric)) %>%
                dplyr::mutate(Gene_Name = stringr::str_to_title(Gene_Name))
            # Grab contrast‑specific columns
            lfc_col <- paste0(contrast, "_diff")
            padj_col <- paste0(contrast, "_p.adj")

            if (!all(c(lfc_col, padj_col) %in% colnames(rd))) {
                stop("contrast columns missing")
            }

            rd <- rd |>
                dplyr::transmute(Gene_ID, Gene_Name, diff = .data[[lfc_col]], p.adj = .data[[padj_col]])
            # message("Extracted features for ", datatype, ": ", nrow(rd), " rows\n COLS: ", paste(colnames(rd), collapse = ", "), "\nHEAD: ", head(rd, n = 3))
            rd
        }

        normalize_ids <- function(dep, datatype) {
            if (methods::is(dep, "DEGdata") || datatype == "rnaseq") {
                rd <- as.data.frame(dep@test_result)
                # If gene IDs/names not separately present, split rownames once
                if (!("Gene_ID" %in% colnames(rd))) {
                    rd$Gene_ID <- sub("^.*_", "", rownames(rd))
                    rd$Gene_Name <- sub("_.*$", "", rownames(rd))
                }
                # message("Normalizing IDs for datatype ", datatype, ": ", nrow(rd), " rows\n COLS: ", paste(colnames(rd), collapse = ", "), "\nHEAD: ", head(rd, n = 1))
            } else {
                rd <- as.data.frame(SummarizedExperiment::rowData(dep))
                # message("Normalizing IDs for datatype ", datatype, ": ", nrow(rd), " rows\n COLS: ", paste(colnames(rd), collapse = ", "), "\nHEAD: ", head(rd, n = 1))
                # Often rowData already has Gene_ID & Gene_Name
                if (!("Gene_ID" %in% colnames(rd))) {
                    stop("No Gene_ID found in rowData for datatype: ", datatype)
                }
            }

            # Ensure correctly typed/clean
            rd <- rd %>%
                dplyr::mutate(
                    Gene_ID   = as.character(Gene_ID),
                    Gene_Name = stringr::str_to_title(as.character(Gene_Name))
                ) %>%
                dplyr::select(where(~ !is.numeric(.)), where(is.numeric))

            rd
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

        get_depflt_matrix <- function(dep_obj, lfc_cut, p_cut, cids) {
            dep_obj <- strip_sig(dep_obj)
            sig <- DEP2::add_rejections(dep_obj, alpha = p_cut, lfc = lfc_cut)
            # If there are no significant rows, or the object is empty, exit early
            if (is.null(sig) || nrow(sig) == 0L) {
                return(list(mat = NULL, mat_scaled = NULL))
            }
            mat <- as.data.frame(SummarizedExperiment::assay(sig))

            if (methods::is(sig, "DEGdata")) {
                df <- as.data.frame(sig@test_result)
                df$id_name <- rownames(df)
                # Rownames are format ENDARG000001_slc01
                # Handle cases where there might not be an underscore
                df$Gene_Name <- str_split_fixed(df$id_name, "_", 2)[,1]
                df$Gene_ID   <- str_split_fixed(df$id_name, "_", 2)[,2]
                df$XID <- "" # There are no XID in rnaseq data
                df <- cbind(df, mat)
            } else {
                df <- cbind(as.data.frame(SummarizedExperiment::rowData(sig)), mat)  
            }      
            mat <- mat[which(df$Gene_ID %in% cids), , drop = FALSE]
            # If assay is NULL or empty, exit early
            mat <- as.matrix(mat)
            df_filt <- df[which(df$Gene_ID %in% cids), , drop = FALSE]
            rownames(mat) <- paste0(df_filt$Gene_ID, "_", df_filt$Gene_Name)
            if (is.null(mat) || length(mat) == 0L) {
                return(list(mat = NULL, mat_scaled = NULL))
            }
            # sanitize matrix (rows/cols need finite variance and ≥2 finite values)
            row_ok <- rowSums(is.finite(mat)) >= 2 & apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
            col_ok <- colSums(is.finite(mat)) >= 2 & apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
            mat <- mat[row_ok, col_ok, drop = FALSE]            
            mat_scaled <- safe_row_scale(mat)
            message("Scaled matrix: ", nrow(mat_scaled), " rows, ", ncol(mat_scaled), " cols", "\n IDS: ", paste(colnames(mat_scaled), collapse = ", "), "\nGENES: ", head(rownames(mat_scaled), n = 3), "\n")
            mat_scaled
        }

        get_df_from_dep <- function(dep) {            
            if (methods::is(dep, "DEGdata")) {
                res <- as.data.frame(dep@test_result) 
                res$id_name <- rownames(res)
                # Rownames are format ENDARG000001_slc01
                res$Gene_Name <- str_split_fixed(res$id_name, "_", 2)[,1]
                res$Gene_ID   <- str_split_fixed(res$id_name, "_", 2)[,2]
                res$XID <- "" # There are no XID in rnaseq data
            } else {
                res <- as.data.frame(SummarizedExperiment::rowData(dep))
            }
            return(res)
        }

        get_matrix_from_dep <- function(dep) {    
            mat <- as.data.frame(SummarizedExperiment::assay(dep))
            if (methods::is(dep, "DEGdata")) {
                df <- as.data.frame(dep@test_result)
                df$id_name <- rownames(df)
                # Rownames are format ENDARG000001_slc01
                df$Gene_Name <- str_split_fixed(df$id_name, "_", 2)[,1]
                df$Gene_ID   <- str_split_fixed(df$id_name, "_", 2)[,2]
                df$XID <- "" # There are no XID in rnaseq data
            } else {
                df <- as.data.frame(SummarizedExperiment::rowData(dep))
            }       
            #mat$names <- rownames(mat)
            mat$Gene_Name <- stringr::str_to_title(df$Gene_Name)  # gsub("_.*", "", mat$names, perl = TRUE)
            mat$Gene_ID <- df$Gene_ID  # gsub(".*_", "", mat$names, perl = TRUE)
            if("pepG" %in% colnames(df)) {
                mat$pepG <- df$pepG
            }
            if("Protein_ID" %in% colnames(df)) {
                mat$Protein_ID <- df$Protein_ID
            }
            mat$XID <- df$XID
            return(mat)
        }

        # Filter each table by thresholds
        for (tbl in selected_tables) {
            contrast <- input[[paste0("pi_comparison_selected_", tbl)]]
            dep <- rv$dep_output[[tbl]]
            datatype <- rv$datatype[[tbl]]

            feats <- extract_features(dep, datatype, contrast)            
            keep <- with(
                feats,
                !is.na(Gene_ID) &
                    !is.na(diff) &
                    !is.na(p.adj) &
                    abs(diff) >= input$lfc_thresh_pi &
                    p.adj <= input$pval_thresh_pi
            )

            filtered_ids[[tbl]] <- unique(feats$Gene_ID[keep])

            dim_info[[tbl]] <- list(
                original   = dim(SummarizedExperiment::assay(dep)),
                filtered   = c(length(filtered_ids[[tbl]]), ncol(SummarizedExperiment::assay(dep)))
            )
        }

        # Intersect IDs
        common_ids <- Reduce(intersect, filtered_ids)
        if (length(common_ids) == 0) {
            showNotification("No intersecting significant IDs found.", type = "error")
            rv$intersected_tables_processed <- NULL
            rv$integration_preview_dims <- NULL
            return()
        }

        # Subset SummarizedExperiment::assays by intersected gene names
        intersected_list <- lapply(selected_tables, function(tbl) {
            dep <- rv$dep_output[[tbl]]

            res <- get_df_from_dep(dep)
            mat <- get_matrix_from_dep(dep)
            dep_flt <- res[res$Gene_ID %in% common_ids | res$XID %in% common_ids, ]  
            mat_flt <- mat[mat$Gene_ID %in% common_ids | mat$XID %in% common_ids, ] 

            data <- dplyr::inner_join(mat_flt, dep_flt, by = "Gene_ID", keep=FALSE, suffix=c("",".y")) %>%
                select(-ends_with(".y"))

            # Generate unique IDs
            if ("pepG" %in% colnames(data)) {
                data$unique_id <- paste(data$Gene_Name, data$pepG, sep = "_")
            } else if ("Protein_ID" %in% colnames(data)) {
                data$unique_id <- paste(data$Gene_Name, data$Protein_ID, sep = "_")
            } else {
                data$unique_id <- data$Gene_Name
            }            
            data
        })

        names(intersected_list) <- selected_tables

        # Extract intersected expression matrices per table
        intersected_matrix <- lapply(selected_tables, function(tbl) {
            dep <- rv$dep_output[[tbl]]
            datatype <- rv$datatype[[tbl]]
            mat <- get_depflt_matrix(dep, input$lfc_thresh_pi, input$pval_thresh_pi, common_ids)
            mat
        })
        names(intersected_matrix) <- selected_tables
        
        #message("Intersected matrices extracted: ", paste(sapply(selected_tables, function(tbl) {
        #    if (is.null(intersected_matrix[[tbl]])) {
        #        paste0(tbl, "NULL")
        #    } else {
        #        paste(tbl, dim(intersected_matrix[[tbl]]), collapse = " x ")
        #    }
        #}), collapse = ", "), "\n")

        # Record dimensions
        for (tbl in selected_tables) {
            dim_info[[tbl]]$intersected <- dim(intersected_matrix[[tbl]])
        }

        # Prepare scatter data: wide table by Gene_ID
        scatter_data <- NULL
        for (tbl in selected_tables) {
            contrast <- input[[paste0("pi_comparison_selected_", tbl)]]
            dep <- rv$dep_output[[tbl]]
            datatype <- rv$datatype[[tbl]]
            feats <- extract_features(dep, datatype, contrast)
            #message("Features for scatter from ", tbl, ": ", nrow(feats), " rows\n COLS: ", paste(colnames(feats), collapse = ", "), "\nHEAD: ", head(feats$Gene_ID, n = 3))
            feats <- feats[feats$Gene_ID %in% common_ids, ]
            feats <- feats |>
                dplyr::select(Gene_ID, Gene_Name, diff) |>
                dplyr::rename(!!paste0("LFC_", tbl) := diff)

            scatter_data <- if (is.null(scatter_data)) feats else dplyr::full_join(scatter_data, feats, by = c("Gene_ID", "Gene_Name"))
        }

        # Build pairwise scatterplots
        lfc_cols <- grep("^LFC_", names(scatter_data), value = TRUE)
        plot_list <- list()
        for (i in 1:(length(lfc_cols) - 1)) {
            for (j in (i + 1):length(lfc_cols)) {
                df_plot <- scatter_data |>
                    dplyr::select(Gene_Name, x = !!lfc_cols[i], y = !!lfc_cols[j]) |>
                    na.omit()

                if (nrow(df_plot) > 0) {
                    p <- ggplot(df_plot, aes(x, y, text = Gene_Name)) +
                        geom_point(alpha = .7, color = "#2b8cbe") +
                        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
                        labs(title = paste(lfc_cols[i], "vs", lfc_cols[j]), x = paste(lfc_cols[i]), y = paste(lfc_cols[j])) +
                        theme_minimal()
                    plot_list[[paste(lfc_cols[i], lfc_cols[j], sep = "_vs_")]] <- p
                }
            }
        }

        message("Processed integration with ", length(common_ids), " common IDs across ", length(selected_tables), " tables.\nIntersect: ", paste(selected_tables, collapse = ", "), "\nDimensions: ", paste(sapply(selected_tables, function(tbl) {
            dim_info[[tbl]]$original[1]
        }), collapse = " x "), " original rows; ", paste(sapply(selected_tables, function(tbl) {
            dim_info[[tbl]]$filtered[1]
        }), collapse = " x "), " filtered rows; ", paste(sapply(selected_tables, function(tbl) {
            dim_info[[tbl]]$intersected[1]
        }), collapse = " x "), " intersected rows.\n")

        # Save results
        rv$scatter_plots <- plot_list
        rv$intersected_tables_processed <- intersected_list
        rv$intersected_matrix_processed <- intersected_matrix
        rv$integration_preview_dims <- lapply(selected_tables, function(tbl) dim_info[[tbl]])
        names(rv$integration_preview_dims) <- selected_tables
        # Estimate cluster number on stacked data
        # Check if all matrices have the same number of columns
        col_counts <- sapply(intersected_matrix, ncol)
        all_same_cols <- length(unique(col_counts)) == 1

        # Conditionally create stacked
        if (all_same_cols) {
            stacked <- do.call(rbind, intersected_matrix)  # genes x samples 
            # If same number of columns stacking is possible
        } else {
            stacked <- intersected_matrix[[1]]  # use only the first table
            # If differing number of columns use only first matrix for optimal k computation
        }
        if (nrow(stacked) < 2) {
            showNotification("Not enough genes for clustering.", type = "error")
            rv$optimal_k <- NULL
            return()
        }
        optimal_k <- safe_nbclust(stacked, k_min = 2, k_max = 10)
        rv$optimal_k <- ifelse(is.na(optimal_k), NULL, optimal_k)
    })
}
