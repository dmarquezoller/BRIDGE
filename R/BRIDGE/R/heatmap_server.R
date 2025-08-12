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
            raw_input <- list(
                raw_data = isolate(rv$tables[[tbl_name]]),
                timepoint_cols = isolate(rv$time_cols[[tbl_name]]), 
                datatype = isolate(rv$datatype[[tbl_name]])
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


#' Server module for DEP heatmap + volcano
#' Expects:
#'   rv$tables[[tbl]], rv$time_cols[[tbl]], rv$datatype[[tbl]]
#'   rv$ht_matrix[[tbl]], rv$optimal_k_individual[[tbl]]
#' @export
dep_heatmap_server <- function(input, output, session, rv, cache) {

  # 0) Register everything once we know the table names
  shiny::observeEvent(rv$table_names, {
    tbls <- isolate(rv$table_names)
    req(length(tbls) > 0)
    message("dep_heatmap_server: registering handlers for tables: ",
            paste(tbls, collapse = ", "))

    lapply(tbls, function(.tbl) local({
      tbl_name <- .tbl
      message("[", tbl_name, "] init")

      # ---- 1) One-time UI setup when data appears ----
      shiny::observeEvent(rv$tables[[tbl_name]], ignoreInit = FALSE, once = TRUE, {
        message("[", tbl_name, "] UI selectize init")
        updateSelectizeInput(
          session, paste0("volcano_search_", tbl_name),
          choices = rv$tables[[tbl_name]]$Gene_Name, server = TRUE
        )
        updateSelectizeInput(
          session, paste0("search_gene_", tbl_name),
          choices = rv$tables[[tbl_name]]$Gene_Name, server = TRUE
        )
      })

      # ---- 2) Base DEP output per table (compute/load) ----
      dep_base <- shiny::reactiveVal(NULL)  # SummarizedExperiment with stats etc.

      # Cache key depends only on data/time cols/datatype (NOT on cutoffs)
      base_key <- shiny::reactive({
        cols_key <- paste(isolate(rv$time_cols[[tbl_name]]), collapse = "_")
        paste("DEPBASE", tbl_name, cols_key, isolate(rv$datatype[[tbl_name]]), sep = "_")
      })

      shiny::observeEvent(
        list(rv$tables[[tbl_name]], rv$time_cols[[tbl_name]], rv$datatype[[tbl_name]]),
        {
          key <- base_key()
          dt  <- isolate(rv$datatype[[tbl_name]])
          dat <- isolate(rv$tables[[tbl_name]])
          req(!is.null(dat), !is.null(dt))

          if (cache$exists(key)) {
            message("[", tbl_name, "] load DEP base from cache: ", key)
            dep_out <- cache$get(key)
          } else {
            message("[", tbl_name, "] compute DEP base: ", key, " (", dt, ")")
            dep_out <- switch(dt,
              proteomics        = dep2_proteomics(dat, tbl_name, rv),
              phosphoproteomics = dep2_phosphoproteomics(dat, tbl_name, rv),
              rnaseq            = dep2_rnaseq(dat, tbl_name, rv),
              stop("[", tbl_name, "] unknown datatype: ", dt)
            )
            cache$set(key, dep_out)
          }

          dep_base(dep_out)                  # make it reactive for this table
          rv$dep_output[[tbl_name]] <- dep_out  # keep your existing field updated

          # valid contrasts once per dep_out
          rd_names <- colnames(SummarizedExperiment::rowData(dep_out))
          sig_cols <- grep("_significant$", rd_names, value = TRUE)
          rv$contrasts[[tbl_name]] <- sub("_significant$", "", sig_cols)

          message("[", tbl_name, "] DEP base ready; contrasts: ",
                  paste(rv$contrasts[[tbl_name]], collapse = ", "))
        },
        ignoreInit = FALSE
      )

      # =====================================================
      # 3) VOLCANO — heavy compute -> future, store to rv
      # =====================================================
      volcano_res <- shiny::reactiveVal(NULL)

      shiny::observeEvent(input[[paste0("compute_volcano_", tbl_name)]], {
        message("[", tbl_name, "] compute_volcano clicked")
        dep_out <- dep_base(); shiny::req(dep_out)  # ensure base is there

        deps <- list(
          contrast   = isolate(input[[paste0("comparison_volcano_", tbl_name)]]),
          pcut       = 10^(-isolate(input[[paste0("volcano_pcutoff_", tbl_name)]])),
          fccut      = isolate(input[[paste0("volcano_fccutoff_", tbl_name)]]),
          dep_output = dep_out,
          datatype   = isolate(rv$datatype[[tbl_name]])
        )
        shiny::req(!is.null(deps$contrast))

        promises::future_promise({
          lfc_col  <- paste0(deps$contrast, "_diff")
          pval_col <- paste0(deps$contrast, "_p.adj")

          res <- SummarizedExperiment::rowData(
            DEP2::add_rejections(deps$dep_output, alpha = deps$pcut, lfc = deps$fccut)
          )

          if (deps$datatype == "proteomics") {
            df <- data.frame(
              gene_names = stringr::str_to_title(res$Gene_Name),
              pval = res[[pval_col]], log2FC = res[[lfc_col]], stringsAsFactors = FALSE
            )
            df_table <- DEP2::get_results(DEP2::get_signicant(deps$dep_output, deps$contrast))
          } else if (deps$datatype == "phosphoproteomics") {
            df <- data.frame(
              peptide = res$pepG,
              pval = res[[pval_col]], log2FC = res[[lfc_col]], stringsAsFactors = FALSE
            )
            df_table <- DEP2::get_results(DEP2::get_signicant(deps$dep_output, deps$contrast))
          } else {
            df <- data.frame(
              gene_ID = rownames(res),
              pval = res[[pval_col]], log2FC = res[[lfc_col]], stringsAsFactors = FALSE
            )
            df_table <- DEP2::get_results(DEP2::get_signicant(deps$dep_output, deps$contrast))
          }

          list(df = df, table = df_table, pcut = deps$pcut, fccut = deps$fccut, datatype = deps$datatype)
        }, seed = TRUE) %...>% (function(value) {
          message("[", tbl_name, "] volcano future resolved (n=", nrow(value$df), ")")
          volcano_res(value)
        }) %...!% (function(e) {
          showNotification(conditionMessage(e), type = "error")
          message("[", tbl_name, "] volcano error: ", conditionMessage(e))
        })
      }, ignoreInit = TRUE)

      output[[paste0("volcano_", tbl_name)]] <- plotly::renderPlotly({
        res <- volcano_res(); shiny::req(res)

        # lightweight UI bits that should re-draw plot without recomputing:
        hl <- input[[paste0("volcano_search_", tbl_name)]] |>
          stringr::str_split(",", simplify = TRUE) |>
          stringr::str_trim() |>
          stringr::str_to_title()

        if (res$datatype == "proteomics") {
          text_col <- res$df$gene_names
          p <- EnhancedVolcano::EnhancedVolcano(
            res$df, lab = res$df$gene_names, selectLab = c("a"),
            x = "log2FC", y = "pval", title = "",
            pCutoff = res$pcut, FCcutoff = res$fccut,
            pointSize = ifelse(text_col %in% hl, 3, 1),
            legendPosition = "none"
          ) + ggplot2::aes(text = gene_names) + ggplot2::labs(color = "Legend")
        } else if (res$datatype == "phosphoproteomics") {
          text_col <- res$df$peptide
          p <- EnhancedVolcano::EnhancedVolcano(
            res$df, lab = res$df$peptide, selectLab = c("a"),
            x = "log2FC", y = "pval", title = "",
            pCutoff = res$pcut, FCcutoff = res$fccut,
            pointSize = ifelse(text_col %in% hl, 3, 1),
            legendPosition = "none"
          ) + ggplot2::aes(text = peptide) + ggplot2::labs(color = "Legend")
        } else {
          text_col <- res$df$gene_ID
          p <- EnhancedVolcano::EnhancedVolcano(
            res$df, lab = res$df$gene_ID, selectLab = c("a"),
            x = "log2FC", y = "pval", title = "",
            pCutoff = res$pcut, FCcutoff = res$fccut,
            pointSize = ifelse(text_col %in% hl, 3, 1),
            legendPosition = "none"
          ) + ggplot2::aes(text = gene_ID) + ggplot2::labs(color = "Legend")
        }

        plotly::ggplotly(p + ggplot2::aes(x = log2FC, y = -log10(pval)), tooltip = "text")
      }) %>% bindEvent(
        input[[paste0("compute_volcano_", tbl_name)]],          # re-compute
        input[[paste0("volcano_search_", tbl_name)]],           # UI-only redraw
        input[[paste0("volcano_pcutoff_", tbl_name)]],          # UI-only redraw
        input[[paste0("volcano_fccutoff_", tbl_name)]],         # UI-only redraw
        ignoreInit = TRUE
      )

      output[[paste0("volcano_sig_table_", tbl_name)]] <- DT::renderDT({
        res <- volcano_res(); shiny::req(res)
        DT::datatable(
          res$table,
          extensions = "Buttons",
          options = list(scrollX = TRUE, pageLength = 10, dom = "Bfrtip",
                         buttons = c('copy','csv','excel','pdf','print'))
        )
      })

      # =====================================================
      # 4) HEATMAP — heavy compute -> future, store to rv
      # =====================================================
      heat_res <- shiny::reactiveVal(NULL)

      shiny::observeEvent(input[[paste0("recompute_heatmap_", tbl_name)]], {
        message("[", tbl_name, "] recompute_heatmap clicked")
        dep_out <- dep_base(); shiny::req(dep_out)
        ht_mat  <- isolate(rv$ht_matrix[[tbl_name]])
        shiny::req(!is.null(ht_mat))

        stored_k <- isolate(rv$optimal_k_individual[[tbl_name]])

        promises::future_promise({ 
          opt_k <- if (is.null(stored_k)) {
            elbow <- NbClust::NbClust(
              ht_mat, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans"
            )
            as.numeric(names(sort(table(elbow$Best.nc[1, ]), decreasing = TRUE)[1]))
          } else stored_k

          dep_pg_sig <- DEP2::get_signicant(dep_out)
          expr       <- SummarizedExperiment::assay(dep_pg_sig)
          gene_info  <- as.data.frame(SummarizedExperiment::rowData(dep_pg_sig))
          df         <- cbind(gene_info, as.data.frame(expr))
          df         <- df[, c(colnames(gene_info), colnames(expr))]
          df         <- stats::na.omit(df)

          list(optimal_k = opt_k, df = df)
        }, seed = TRUE) %...>% (function(value) {
          message("[", tbl_name, "] heatmap future resolved (rows=", nrow(value$df), ")")
          heat_res(value)
        }) %...!% (function(e) {
          showNotification(conditionMessage(e), type = "error")
          message("[", tbl_name, "] heatmap error: ", conditionMessage(e))
        })
      }, ignoreInit = TRUE)

      output[[paste0("ht_", tbl_name)]] <- shiny::renderPlot({
        dep_out <- dep_base(); shiny::req(dep_out)
        res <- heat_res(); shiny::req(res)

        clustering_enabled <- input[[paste0("clustering_", tbl_name)]]
        num_clusters       <- input[[paste0("num_clusters_", tbl_name)]]

        if (isTRUE(clustering_enabled)) {
          DEP2::plot_heatmap(dep_out, kmeans = TRUE, k = num_clusters)
        } else {
          DEP2::plot_heatmap(dep_out)
        }
      }) %>% bindEvent(
        input[[paste0("recompute_heatmap_", tbl_name)]],        # re-compute
        input[[paste0("clustering_", tbl_name)]],               # UI-only redraw
        input[[paste0("num_clusters_", tbl_name)]],             # UI-only redraw
        ignoreInit = TRUE
      )

      output[[paste0("ht_sig", tbl_name)]] <- DT::renderDT({
        res <- heat_res(); shiny::req(res)
        DT::datatable(
          res$df,
          extensions = "Buttons",
          options = list(scrollX = TRUE, pageLength = 10, dom = "Bfrtip",
                         buttons = c('copy','csv','excel','pdf','print'))
        )
      })

      output[[paste0("optimal_k", tbl_name)]] <- renderUI({
        res <- heat_res(); shiny::req(res)
        if (is.null(rv$optimal_k_individual[[tbl_name]])) {
          rv$optimal_k_individual[[tbl_name]] <- res$optimal_k
        }
        shiny::tags$ul(shiny::tags$li(sprintf(
          "The optimal k for this table following the elbow rule is: %s", res$optimal_k
        )))
      })

    }))  # end local(.tbl)
  }, ignoreInit = FALSE, once = TRUE)  # register once when table_names is known
}