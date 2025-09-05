#' @export
pcaServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        get_gene_loadings <- function(dep, ntop = 500, center = TRUE, scale. = FALSE, pcs = 1:2, rank_by = c("abs", "pos", "neg")) {
            rank_by <- match.arg(rank_by)
            mat <- assay(dep)
            rv <- matrixStats::rowVars(mat)
            select <- head(order(rv, decreasing = TRUE), min(ntop, length(rv)))
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
            pca <- prcomp(t(mat[select, , drop = FALSE]), center = center, scale. = scale.)

            loadings <- pca$rotation # rows = genes (selected); cols = PCs
            var_expl <- pca$sdev^2 / sum(pca$sdev^2)

            df <- as.data.frame(loadings) |>
                rownames_to_column("gene") |>
                pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") |>
                mutate(
                    abs_loading = abs(loading),
                    pc_idx = as.integer(sub("PC", "", PC)),
                    var_explained = var_expl[pc_idx]
                )

            rank_col <- switch(rank_by,
                abs = "abs_loading",
                pos = "loading",
                neg = "loading"
            )
            if (rank_by == "neg") df <- df |> mutate(loading = -loading) # so largest = most negative
            df <- df |>
                group_by(PC) |>
                arrange(desc(!!sym(if (rank_by == "abs") "abs_loading" else "loading"))) |>
                mutate(rank = row_number()) |>
                ungroup()

            # restore original sign if we flipped
            if (rank_by == "neg") df <- df |> mutate(loading = -loading)

            df
        }

        pca_task <- ExtendedTask$new(function(dep) {
            promises::future_promise({
                DESeq2::DESeqTransform(dep)
            })
        })

        pca_loadings_task <- ExtendedTask$new(function(dep) {
            promises::future_promise({
                loadings_df <- get_gene_loadings(dep, ntop = 500, pcs = 1:3, rank_by = "abs")

                contrib <- loadings_df %>%
                    group_by(PC) %>%
                    mutate(contribution = (loading^2) / sum(loading^2)) %>%
                    ungroup()

                # Top 50 genes (by |loading|) for PC1:
                # top_pc1 <- loadings_df |>
                #    filter(PC == "PC1") |>
                #    slice_max(abs_loading, n = 10)
                #
                # top_pc2 <- loadings_df |>
                #    filter(PC == "PC2") |>
                #    slice_max(abs_loading, n = 10)

                top_contrib <- contrib %>%
                    filter(PC == "PC1" | PC == "PC2") %>%
                    arrange(desc(contribution))
                top_contrib
            })
        })

        observeEvent(input$compute,
            {
                req(rv$dep_output[[tbl_name]])
                pca_task$invoke(isolate(rv$dep_output[[tbl_name]]))
                pca_loadings_task$invoke(isolate(rv$dep_output[[tbl_name]]))
            },
            ignoreInit = TRUE
        )

        output$plot <- renderPlot({
            res <- pca_task$result()
            req(res)
            DESeq2::plotPCA(res) +
                ggplot2::ggtitle(paste("PCA for", tbl_name)) +
                ggplot2::theme_minimal()
        })
        output$pcs_panel <- renderUI({
            DT::DTOutput(session$ns("pcs"), height = "300px")
        })

        output$pcs <- DT::renderDT({
            top_contrib <- pca_loadings_task$result()
            if (is.null(top_contrib)) {
                validate(need(FALSE, htmltools::tagList(
                    tags$span(class = "lds-ring", tags$div(), tags$div(), tags$div(), tags$div()),
                    " Loading PCA loadings…"
                )))
            }
            req(!is.null(top_contrib))
            top_contrib <- as.data.frame(top_contrib)
            DT::datatable(top_contrib,
                extensions = "Buttons",
                options = list(
                    scrollX = TRUE, processing = TRUE, pageLength = 10, dom = "Bfrtip",
                    buttons = c("copy", "csv", "excel", "pdf", "print")
                )
            )
        })
    })
}
