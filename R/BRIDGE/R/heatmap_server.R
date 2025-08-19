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
        raw_data       <- raw_input$raw_data
        timepoint_cols <- raw_input$timepoint_cols
        datatype       <- raw_input$datatype

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

    observeEvent(input$compute, {
      ready(TRUE)
      task$invoke(list(
        raw_data       = isolate(rv$tables[[tbl_name]]),
        timepoint_cols = isolate(rv$time_cols[[tbl_name]]),
        datatype       = isolate(rv$datatype[[tbl_name]])
      ))
    }, ignoreInit = TRUE)

    # Mount spinner only after "Compute" is pressed
    output$plot_slot <- renderUI({
      if (!ready()) {
        return(div(style="padding:10px; color:#777;",
                   "Click “Compute Raw Heatmap” to generate the plot."))
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
                              show_row_names = FALSE)
    })
  })
}


#### THIS FUNCTION CONTAINS THE SERVER FOR THE DEP HEATMAP
#' @export
DepHeatmapServer <- function(id, rv, cache, tbl_name) {
  moduleServer(id, function(input, output, session) {

    heatmap_ready <- reactiveVal(FALSE)

    heatmap_task <- ExtendedTask$new(function(ht_inps) {
        promises::future_promise({
            #stored_k   <- ht_inps$stored_k
            dep_output <- ht_inps$dep_output
            mat <- as.matrix(ht_inps$ht_matrix)
            # drop rows with <2 finite values
            ok <- rowSums(is.finite(mat)) >= 2
            mat <- mat[ok, , drop = FALSE]
            
            # open a throwaway graphics device to avoid device warnings (see 3)
            tmp <- tempfile(fileext = ".pdf")
            pdf(tmp); on.exit({ dev.off(); unlink(tmp) }, add = TRUE)

            # Guard against very small matrices
            if (nrow(mat) < 2 || ncol(mat) < 2) {
                return(list(optimal_k = NA_integer_, df = data.frame()))
            }       

            elbow <- NbClust::NbClust(
                mat, distance = "euclidean",
                min.nc = 2, max.nc = 10,
                method = "kmeans",
                index  = "ch"       # <- pick one index; avoids Beale/pf() path
            )
            # Compute optimal_k robustly regardless of Best.nc shape
            best <- elbow$Best.nc
            if (is.matrix(best)) {
                k_vec <- as.integer(best[1, ])
                # mode of candidate ks
                optimal_k <- as.integer(names(which.max(table(k_vec))))
            } else {
                # single index (e.g. "ch") returns a vector
                optimal_k <- as.integer(best[1])
            }

            dep_pg_sig <- DEP2::get_signicant(dep_output)
            expr       <- SummarizedExperiment::assay(dep_pg_sig)
            gene_info  <- as.data.frame(SummarizedExperiment::rowData(dep_pg_sig))
            df         <- cbind(gene_info, as.data.frame(expr))
            df         <- df[, c(colnames(gene_info), colnames(expr)), drop = FALSE]
            df         <- stats::na.omit(df)

            list(optimal_k=optimal_k, df=df)
      }, seed = TRUE)
    })

    observe({
        res <- try(heatmap_task$result(), silent = TRUE)
        if (inherits(res, "try-error")) {
            showNotification(paste("DEP heatmap error:", conditionMessage(attr(res, "condition"))),
                            type = "error")
        }
    })

    observeEvent(input$recompute_heatmap, {
      req(rv$dep_output[[tbl_name]], rv$ht_matrix[[tbl_name]])
      heatmap_ready(TRUE)

      clustering_enabled <- isolate(input$clustering)
      num_clusters       <- isolate(input$num_clusters)
      dep_output         <- isolate(rv$dep_output[[tbl_name]])
      ht_matrix          <- isolate(rv$ht_matrix[[tbl_name]])
      #stored_k           <- isolate(rv$optimal_k_individual[[tbl_name]])
      columns_key        <- paste(isolate(rv$time_cols[[tbl_name]]), collapse = "_")

      key <- if (isTRUE(clustering_enabled)) {
        paste(tbl_name, columns_key, paste0("clusters:", num_clusters), "dep_heatmap", sep = "_")
      } else {
        paste(tbl_name, columns_key, "dep_heatmap", sep = "_")
      }
      rv$current_dep_heatmap_key[[tbl_name]] <- key

      if (!cache$exists(key)) {
        heatmap_task$invoke(list(dep_output=dep_output, ht_matrix=ht_matrix))  #, stored_k=stored_k))
      }
    }, ignoreInit = TRUE)

    get_dep_result <- function() {
      key <- rv$current_dep_heatmap_key[[tbl_name]]
      if (!is.null(key) && cache$exists(key)) return(cache$get(key))
      heatmap_task$result()
    }

    output$ht_slot <- renderUI({
      if (!heatmap_ready()) {
        return(div(style="padding:10px; color:#777;", "Choose settings and click “Recompute heatmap”."))
      }
      shinycssloaders::withSpinner(
        plotOutput(session$ns("ht"), height = "520px"),
        type = 8, color = "#2b8cbe", caption = "Loading..."
      )
    })

    output$ht <- renderPlot({
      req(rv$dep_output[[tbl_name]])
      clustering_enabled <- input$clustering
      num_clusters       <- input$num_clusters
      p_cut   <- input$heatmap_pcutoff
      lfc_cut <- input$heatmap_fccutoff

      res <- get_dep_result(); req(res)
      key <- rv$current_dep_heatmap_key[[tbl_name]]
      if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

      if(is.na(res$optimal_k)){showNotification("The chosen cutoffs leave too few rows in the matrix for sensible optimal k calculation, please consider relaxing your cutoffs")}

      dep_output <- rv$dep_output[[tbl_name]]
      rd <- SummarizedExperiment::rowData(dep_output)
      rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
      SummarizedExperiment::rowData(dep_output) <- rd

      dep_output_filtered <- DEP2::add_rejections(dep_output, alpha = p_cut, lfc = lfc_cut)

      if (isTRUE(clustering_enabled)) {
        DEP2::plot_heatmap(dep_output_filtered, kmeans = TRUE, k = num_clusters)
      } else {
        DEP2::plot_heatmap(dep_output_filtered)
      }    
    })

    output$ht_sig <- DT::renderDT({
      res <- get_dep_result(); req(res)
      key <- rv$current_dep_heatmap_key[[tbl_name]]
      if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

      dep_output   <- rv$dep_output[[tbl_name]]
      p_cut <- input$heatmap_pcutoff; lfc_cut <- input$heatmap_fccutoff

      rd <- SummarizedExperiment::rowData(dep_output)
      rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
      SummarizedExperiment::rowData(dep_output) <- rd    

      #row_order_list <- ComplexHeatmap::row_order(ht)
      gene_names <- SummarizedExperiment::rowData(dep_output)$Gene_Name

      cluster_assignments <- rep(NA, length(gene_names))    
      for (cl in seq_along(row_order_list)) {
        cluster_assignments[row_order_list[[cl]]] <- cl
      }
      res$df$Clusters <- cluster_assignments
      df_filtered <- res$df[res$df$Gene_Name %in% gene_names, , drop = FALSE]

      DT::datatable(df_filtered, extensions="Buttons",
                    options=list(scrollX=TRUE, pageLength=10,
                                 dom="Bfrtip", buttons=c('copy','csv','excel','pdf','print')))
    })

    output$optimal_k <- renderUI({
      res <- get_dep_result(); req(res)
      key <- rv$current_dep_heatmap_key[[tbl_name]]
      if (!is.null(key) && !cache$exists(key)) cache$set(key, res)

      if (is.null(rv$optimal_k_individual[[tbl_name]])) {
        rv$optimal_k_individual[[tbl_name]] <- res$optimal_k
      }
      tags$ul(tags$li(paste("The optimal k for this table (elbow rule):", res$optimal_k)))
    })
  })
}

