#' @export
datalineServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        # --- gene choices for this table ---
        observe({
            req(rv$tables[[tbl_name]])
            # Prefer Gene_Name column; otherwise fall back to rownames
            choices <- rv$tables[[tbl_name]][["Gene_Name"]]
            if (is.null(choices)) choices <- rownames(rv$tables[[tbl_name]])
            updateSelectizeInput(session, "search_gene", choices = sort(unique(choices)), server = TRUE)
        })

        # --- normalization per table/scale ---
        processed_data <- reactive({
            req(rv$tables[[tbl_name]], input$scale)
            raw_data <- rv$tables[[tbl_name]]
            switch(input$scale,
                "Total Intensity" = total_intensity_normalization(raw_data),
                "Continous" = raw_data,
                "Log-scale" = log2_transform(raw_data),
                "Median Normalization" = median_normalization(raw_data),
                "FPKM" = {
                    req(rv$annotation[[tbl_name]])
                    fpkm_normalization(raw_data, rv$annotation[[tbl_name]])
                },
                "TMM" = tmm_normalization(raw_data),
                "CPM" = cpm_normalization(raw_data),
                "TPM" = {
                    req(rv$annotation[[tbl_name]])
                    tpm_normalization(fpkm_normalization(raw_data, rv$annotation[[tbl_name]]))
                },
                raw_data
            )
        })

        genes_chosen <- reactive({
            req(input$search_gene)
            stringr::str_to_lower(trimws(input$search_gene))
        })

        unique_datapoints <- reactive({
            req(rv$data_cols[[tbl_name]])
            rv$data_cols[[tbl_name]] |>
                gsub("_[0-9]+$", "", x = _) |>
                gsub("[_.-]+$", "", x = _) |>
                unique()
        })

        # --- long format for plotting (handles phospho vs others) ---
        long_data <- reactive({
            req(processed_data(), rv$data_cols[[tbl_name]])
            dt <- processed_data()
            tp <- rv$data_cols[[tbl_name]]
            if (identical(rv$datatype[[tbl_name]], "phosphoproteomics")) {
                df <- subset(dt, stringr::str_to_lower(Gene_Name) %in% genes_chosen())
                validate(need(nrow(df) > 0, "No data found for the entered gene(s)."))
                df$Gene_pepG <- paste(stringr::str_to_title(df$Gene_Name), df$pepG, sep = "_")
                keep <- c("Gene_pepG", tp)
                df <- df[, keep, drop = FALSE]
                dl <- tidyr::pivot_longer(df, cols = -Gene_pepG, names_to = "Stage", values_to = "Expression")
                dl$StageGroup <- factor(gsub("[_.-]+$", "", gsub("[0-9]+$", "", dl$Stage)),
                    levels <- unique_datapoints()
                )
                dl
            } else {
                df <- subset(dt, stringr::str_to_lower(Gene_Name) %in% genes_chosen())
                validate(need(nrow(df) > 0, "No data found for the entered gene(s)."))
                keep <- c("Gene_Name", tp)
                df <- df[, keep, drop = FALSE]
                dl <- tidyr::pivot_longer(df, cols = -Gene_Name, names_to = "Stage", values_to = "Expression")
                dl$StageGroup <- factor(gsub("[_.-]+$", "", gsub("[0-9]+$", "", dl$Stage)),
                    levels <- unique_datapoints()
                )
                dl
            }
        })

        # --- table output (wide) ---
        output$data_plot_dt <- DT::renderDT({
            df <- processed_data()
            tp <- rv$data_cols[[tbl_name]]
            if (identical(rv$datatype[[tbl_name]], "phosphoproteomics")) {
                keep <- c("Gene_Name", "pepG", tp)
            } else {
                keep <- c("Gene_Name", tp)
            }
            DT::datatable(df[, keep, drop = FALSE],
                extensions = "Buttons",
                options = list(
                    scrollX = TRUE, pageLength = 10,
                    dom = "Bfrtip",
                    buttons = c("copy", "csv", "excel", "pdf", "print")
                )
            )
        })

        # --- plot output ---
        output$data_plot <- renderPlot({
            dl <- long_data()
            if (identical(rv$datatype[[tbl_name]], "phosphoproteomics")) {
                data_avg <- dplyr::summarise(dplyr::group_by(dl, StageGroup, Gene_pepG),
                    MeanExpression = mean(Expression, na.rm = TRUE),
                    .groups = "drop"
                )
                p <- ggplot2::ggplot() +
                    ggplot2::geom_point(
                        data = dl,
                        ggplot2::aes(
                            x = StageGroup, y = Expression,
                            color = as.factor(stringr::str_to_title(Gene_pepG))
                        ),
                        size = 3, alpha = 0.6
                    ) +
                    ggplot2::geom_line(
                        data = data_avg,
                        ggplot2::aes(
                            x = StageGroup, y = MeanExpression,
                            color = as.factor(stringr::str_to_title(Gene_pepG)),
                            group = stringr::str_to_title(Gene_pepG)
                        ),
                        linewidth = 1.2
                    ) +
                    ggplot2::geom_point(
                        data = data_avg,
                        ggplot2::aes(
                            x = StageGroup, y = MeanExpression,
                            color = as.factor(stringr::str_to_title(Gene_pepG))
                        ),
                        size = 4, shape = 17
                    ) +
                    ggplot2::labs(
                        x = "Stage",
                        y = "log Expression",
                        title = paste(
                            "Peptide Expression data Line",
                            stringr::str_to_title(input$search_gene)
                        ),
                        color = "Proteins"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(
                        text = ggplot2::element_text(size = 14),
                        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                    )
                if (identical(input$scale, "Log-scale")) p <- p + ggplot2::scale_y_log10()
                p
            } else {
                data_avg <- dplyr::summarise(dplyr::group_by(dl, StageGroup, Gene_Name),
                    MeanExpression = mean(Expression, na.rm = TRUE),
                    .groups = "drop"
                )
                ggplot2::ggplot() +
                    ggplot2::geom_point(
                        data = dl,
                        ggplot2::aes(
                            x = StageGroup, y = Expression,
                            color = as.factor(stringr::str_to_title(Gene_Name))
                        ),
                        size = 3, alpha = 0.6
                    ) +
                    ggplot2::geom_line(
                        data = data_avg,
                        ggplot2::aes(
                            x = StageGroup, y = MeanExpression,
                            color = as.factor(stringr::str_to_title(Gene_Name)),
                            group = stringr::str_to_title(Gene_Name)
                        ),
                        linewidth = 1.2
                    ) +
                    ggplot2::geom_point(
                        data = data_avg,
                        ggplot2::aes(
                            x = StageGroup, y = MeanExpression,
                            color = as.factor(stringr::str_to_title(Gene_Name))
                        ),
                        size = 4, shape = 17
                    ) +
                    ggplot2::labs(
                        x = "Stage",
                        y = "log Expression",
                        title = paste(
                            "Protein Expression data Line",
                            stringr::str_to_title(input$search_gene)
                        ),
                        color = "Proteins"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(
                        text = ggplot2::element_text(size = 14),
                        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                    )
            }
        })
    })
}
