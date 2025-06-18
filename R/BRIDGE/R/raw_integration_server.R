#' @export
raw_integration <- function(input, output, session, rv, combined_data){

    shiny::observeEvent(input$integrate_data, {
        req(length(input$integration) > 0)
        selected_tables <- input$integration
        selected_types <- sapply(selected_tables, function(tbl) rv$datatype[[tbl]])

        # Extract selected columns
        selected_lists <- lapply(selected_tables, function(tbl) input[[paste0("cols_selected_int", tbl)]])
        
        # Validate selection
        if (any(sapply(selected_lists, is.null))) {
        combined_data(data.frame(Message = "Select columns from each table before integrating."))
        return()
        }
        if (length(unique(sapply(selected_lists, length))) > 1) {
        combined_data(data.frame(Message = "Error: Column numbers differ. Matching not possible."))
        return()
        }

        # Determine if phospho integration pipeline is needed
        use_phospho_pipeline <- "phosphoproteomics" %in% selected_types

        reference_time_names <- selected_lists[[1]]
        tables_to_join <- list()

        for (i in seq_along(selected_tables)) {
        tbl <- selected_tables[[i]]
        selected_time_cols <- selected_lists[[i]]
        df <- rv$tables[[tbl]][, c(rv$id_cols[[tbl]], selected_time_cols), drop = FALSE]

        # Determine if phospho
        is_phospho <- rv$datatype[[tbl]] == "phosphoproteomics"

        # Construct ID column
        if (is_phospho && all(c("Gene_Name", "pepG") %in% names(df))) {
            df$unique_id <- paste(df$Gene_Name, df$pepG, sep = "_")
        } else {
            df$unique_id <- df$Gene_Name
        }

        # Ensure ID columns match across tables
        id_cols <- setdiff(c("unique_id", rv$id_cols[[tbl]]), selected_time_cols)
        all_id_names <- unique(unlist(lapply(selected_tables, function(x) c("unique_id", rv$id_cols[[x]]))))

        # Add missing ID columns with NA
        missing_ids <- setdiff(all_id_names, names(df))
        for (col in missing_ids) df[[col]] <- NA

        # Rename time columns (to match the first table)
        if (i > 1) {
            colnames(df)[match(selected_time_cols, names(df))] <- reference_time_names
        }

        # Reorder columns
        df <- df[, c(all_id_names, reference_time_names), drop = FALSE]
        df$source <- tbl
        tables_to_join[[tbl]] <- df
        }

        # Combine all
        combined_df <- dplyr::bind_rows(tables_to_join)
        combined_data(combined_df)

        # Update search box
        updateSelectizeInput(session, "search_gene_integration", choices = unique(combined_df$unique_id), server = TRUE)

        # Render multi-timeline plot
        int_timeline_server(input, output, session, combined_df, reference_time_names)
    })
}