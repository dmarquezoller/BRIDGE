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
raw_heatmap_server <- function(input, output, session, rv) {
    shiny::observe({
    lapply(rv$table_names, function(tbl_name) {

        raw_heatmap_plot <- ExtendedTask$new(function(raw_input){
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

        shiny::observeEvent(input[[paste0("compute_raw_ht_", tbl_name)]], {
  
          raw_data <- isolate(rv$tables[[tbl_name]])
          timepoint_cols <- isolate(rv$time_cols[[tbl_name]])
          datatype <- isolate(rv$datatype[[tbl_name]])

          raw_input <- list(
            raw_data = raw_data,
            timepoint_cols = timepoint_cols, 
            datatype = datatype
          )

          raw_heatmap_plot$invoke(raw_input)
        })
          
        output[[paste0("raw_ht_", tbl_name)]] <- shiny::renderPlot({
            res <- raw_heatmap_plot$result()
            shiny::req(res)
            ComplexHeatmap::Heatmap(res$clean_matrix, show_row_dend = FALSE, show_row_names = FALSE)
        })       
    })
  })
}


#### THIS FUNCTION CONTAINS BOTH THE SERVER FOR THE DEP HEATMAP AND THE VOLCANO (because both use the dep_output)
#' @export
dep_heatmap_server <- function(input, output, session, rv, cache) {
  shiny::observe({
    lapply(rv$table_names, function(tbl_name) {
      updateSelectizeInput(session, paste0("volcano_search_",tbl_name), choices = rv$tables[[tbl_name]]$Gene_Name, server = TRUE)
      updateSelectizeInput(session, paste0("search_gene_", tbl_name), choices = rv$tables[[tbl_name]]$Gene_Name, server = TRUE)
      time_cols_hash <- rv$time_cols[[tbl_name]]
      cache_key <- paste(tbl_name, paste(rv$time_cols[[tbl_name]], collapse = "_"), "dep", sep = "_")
      print(paste("Running analysis for:", tbl_name, "datatype:", rv$datatype[[tbl_name]]))

      # Recompute heatmap if the button is clicked

      if (cache$exists(cache_key)) {
        message("Loading DEP output from cache: ", cache_key)
        dep_output <- cache$get(cache_key)
      } else {
        message("Computing and caching DEP output: ", cache_key)
        if (rv$datatype[[tbl_name]] == 'proteomics') {
          dep_output <- dep2_proteomics(rv$tables[[tbl_name]], tbl_name, rv)
        } else if (rv$datatype[[tbl_name]] == 'phosphoproteomics') {
          dep_output <- dep2_phosphoproteomics(rv$tables[[tbl_name]], tbl_name, rv)
        } else if (rv$datatype[[tbl_name]] == 'rnaseq') {
          dep_output <- dep2_rnaseq(rv$tables[[tbl_name]], tbl_name, rv)
        }
        cache$set(cache_key, dep_output)
      }

     shiny::req(input[[paste0("volcano_pcutoff_", tbl_name)]],  input[[paste0("volcano_fccutoff_", tbl_name)]])
      pcut <- input[[paste0("volcano_pcutoff_", tbl_name)]]
      lfcut <- input[[paste0("volcano_fccutoff_", tbl_name)]]

      dep_output <- DEP2::add_rejections(dep_output, alpha = 10^-pcut, lfc = lfcut)
      rd_names <- colnames( SummarizedExperiment::rowData(dep_output))
      sig_cols <- grep("_significant$", rd_names, value = TRUE)
      valid_contrasts <- sub("_significant$", "", sig_cols)  
      #Isolate function to avoid retriggering the shiny::observe block
      isolate({ 
        rv$dep_output[[tbl_name]] <- dep_output
        rv$contrasts[[tbl_name]] <- valid_contrasts
      })  

#### TESTING NON BLOCKING

    volcano_plot <- ExtendedTask$new(function(deps) {
        promises::future_promise({
            dep_output <- deps$dep_output
            contrast   <- deps$contrast
            datatype   <- deps$datatype
            pcut       <- deps$pcut
            fccut      <- deps$fccut

            lfc_col <- paste0(contrast, "_diff")
            pval_col <- paste0(contrast, "_p.adj")
            res <- SummarizedExperiment::rowData(dep_output)

            if (datatype == "proteomics") {
            df <- data.frame(
                gene_names = stringr::str_to_title(res$Gene_Name),
                pval       = res[[pval_col]],
                log2FC     = res[[lfc_col]],
                stringsAsFactors = FALSE
            )
            sig_se   <- DEP2::get_signicant(dep_output, contrast)
            df_table <- DEP2::get_results(sig_se)
            } else if (datatype == "phosphoproteomics") {
            df <- data.frame(
                peptide = res$pepG,
                pval    = res[[pval_col]],
                log2FC  = res[[lfc_col]],
                stringsAsFactors = FALSE
            )
            sig_se   <- DEP2::get_signicant(dep_output, contrast)
            df_table <- DEP2::get_results(sig_se)
            } else { # rnaseq
            df <- data.frame(
                gene_ID = rownames(res),
                pval    = res[[pval_col]],
                log2FC  = res[[lfc_col]],
                stringsAsFactors = FALSE
            )
            df_table <- DEP2::get_results(DEP2::get_signicant(dep_output, contrast))
            }

            list(df = df, table = df_table, pcut = pcut, fccut = fccut, datatype = datatype)
        })
    })

    shiny::observeEvent(input[[paste0("compute_volcano_", tbl_name)]], {
        contrast <- isolate(input[[paste0("comparison_volcano_", tbl_name)]])
        pcut <- isolate(10^(-input[[paste0("volcano_pcutoff_", tbl_name)]]))
        fccut <- isolate(input[[paste0("volcano_fccutoff_", tbl_name)]])
        highlight <- isolate(input[[paste0("volcano_search_", tbl_name)]]) %>%
          str_split(",", simplify = TRUE) %>%
          str_trim() %>%
          stringr::str_to_title()
        dep_output <- isolate(rv$dep_output[[tbl_name]])
        datatype <- isolate(rv$datatype[[tbl_name]])

        deps <- list(
          contrast = contrast,
          pcut = pcut,
          fccut = fccut,
          highlight = highlight,
          dep_output = dep_output,
          datatype = datatype
        )

        volcano_plot$invoke(deps)
    })

    output[[paste0("volcano_", tbl_name)]] <- plotly::renderPlotly({
        res <- volcano_plot$result()
        shiny::req(res)

        if (res$datatype == "proteomics") {
            p <- EnhancedVolcano::EnhancedVolcano(
            res$df, lab = res$df$gene_names, x = "log2FC", y = "pval",
            title = "", pCutoff = res$pcut, FCcutoff = res$fccut, legendPosition = "none"
            )
        } else if (res$datatype == "phosphoproteomics") {
            p <- EnhancedVolcano::EnhancedVolcano(
            res$df, lab = res$df$peptide, x = "log2FC", y = "pval",
            title = "", pCutoff = res$pcut, FCcutoff = res$fccut, legendPosition = "none"
            )
        } else {
            p <- EnhancedVolcano::EnhancedVolcano(
            res$df, lab = res$df$gene_ID, x = "log2FC", y = "pval",
            title = "", pCutoff = res$pcut, FCcutoff = res$fccut, legendPosition = "none"
            )
        }

        # Build plotly on main
        plotly::ggplotly(p + ggplot2::aes(x = log2FC, y = -log10(pval)), tooltip = "text")
    })

    output[[paste0("volcano_sig_table_", tbl_name)]] <- DT::renderDT({
        res <- volcano_plot$result()
        shiny::req(res)
        DT::datatable(res$table, options = list(scrollX = TRUE, pageLength = 10))
    })



#### END OF TESTING NON BLOCKING

    # One reactiveVal per table to hold the latest result
    dep_res <- reactiveVal(NULL)

    plot_dep_heatmap <- ExtendedTask$new(function(ht_inps) {
        promises::future_promise({
            ht_matrix  <- ht_inps$ht_matrix
            dep_output <- ht_inps$dep_output
            stored_k   <- ht_inps$stored_k
            optimal_k <- if (is.null(stored_k)) {
            elbow <- NbClust::NbClust(ht_matrix, distance = "euclidean",
                                        min.nc = 2, max.nc = 10, method = "kmeans")
            as.numeric(names(sort(table(elbow$Best.nc[1, ]), decreasing = TRUE)[1]))
            } else stored_k

            dep_pg_sig <- DEP2::get_signicant(dep_output)
            expr       <- SummarizedExperiment::assay(dep_pg_sig)
            gene_info  <- as.data.frame(SummarizedExperiment::rowData(dep_pg_sig))
            df         <- cbind(gene_info, as.data.frame(expr))
            df         <- df[, c(colnames(gene_info), colnames(expr))]

            list(optimal_k = optimal_k, df = df)
        })
    })


    observeEvent(input[[paste0("recompute_heatmap_", tbl_name)]], {
        clustering_enabled <- isolate(input[[paste0("clustering_", tbl_name)]])
        num_clusters       <- isolate(input[[paste0("num_clusters_", tbl_name)]])
        dep_output         <- isolate(rv$dep_output[[tbl_name]])
        ht_matrix          <- isolate(rv$ht_matrix[[tbl_name]])
        stored_k           <- isolate(rv$optimal_k_individual[[tbl_name]])
        columns_key        <- paste(rv$time_cols[[tbl_name]], collapse = "_")

        key <- if (isTRUE(clustering_enabled)) {
            paste(tbl_name, columns_key, paste0("clusters:", num_clusters), "dep_heatmap", sep = "_")
        } else {
            paste(tbl_name, columns_key, "dep_heatmap", sep = "_")
        }
        rv$current_dep_heatmap_key[[tbl_name]] <- key  # stash the key for renderers

        if (cache$exists(key)) {
            # Cached: nothing to compute; renderers will read from cache below.
            message("Loading DEP heatmap from cache: ", key)
        } else {
            message("Computing and caching DEP heatmap: ", key)
            plot_dep_heatmap$invoke(list(
            dep_output = dep_output,
            ht_matrix  = ht_matrix,
            stored_k   = stored_k
            ))
        }
    }, ignoreInit = TRUE)

    # Helper to fetch the value either from cache or from the task result
    get_dep_result <- function(tbl_name) {
        key <- rv$current_dep_heatmap_key[[tbl_name]]
        if (!is.null(key) && cache$exists(key)) {
            return(cache$get(key))
        }
        # not cached yet — if the task finished, result() will be non-NULL
        plot_dep_heatmap$result()
    }

    # Renderers: read result reactively; write to cache on the main thread once
    output[[paste0("ht_", tbl_name)]] <- shiny::renderPlot({
        # Ensure inputs used to build the plot are present (and for cache key)
        clustering_enabled <- input[[paste0("clustering_", tbl_name)]]
        num_clusters       <- input[[paste0("num_clusters_", tbl_name)]]
        dep_output         <- rv$dep_output[[tbl_name]]
        req(dep_output)  # avoid early NULLs

        res <- get_dep_result(tbl_name)
        req(res)  # wait until we actually have data

        # Cache-once guard
        key <- rv$current_dep_heatmap_key[[tbl_name]]
        if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

        # Build plot on main process
        if (isTRUE(clustering_enabled)) {
            DEP2::plot_heatmap(dep_output, kmeans = TRUE, k = num_clusters)
        } else {
            DEP2::plot_heatmap(dep_output)
        }
    })

    output[[paste0("ht_sig", tbl_name)]] <- DT::renderDT({
        res <- get_dep_result(tbl_name)
        req(res)
        key <- rv$current_dep_heatmap_key[[tbl_name]]
        if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

        DT::datatable(res$df, options = list(scrollX = TRUE, pageLength = 10))
    })

    output[[paste0("optimal_k", tbl_name)]] <- renderUI({
        res <- get_dep_result(tbl_name)
        req(res)
        key <- rv$current_dep_heatmap_key[[tbl_name]]
        if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

        if (is.null(rv$optimal_k_individual[[tbl_name]])) {
            rv$optimal_k_individual[[tbl_name]] <- res$optimal_k
        }
        htmltools::tagList(
            shiny::tags$ul(
            shiny::tags$li(paste("The optimal k for this table following the elbow rule is:", res$optimal_k))
            )
        )
    }) # End of lapply for each table
  })
  }) # End of moduleServer
}