source("R/global.R")

pdf(file = NULL)

heatmap_server <- function(id, rv) {
  
  moduleServer(id, function(input, output, session) {
    observe({
      lapply(rv$table_names, function(tbl_name) {
        output[[paste0("raw_ht_", tbl_name)]] <- renderPlot({
          Heatmap(rv$ht_matrix[[tbl_name]], show_row_dend = FALSE, show_row_names = FALSE)
        })
      })
    })
  })
}

raw_heatmap_server <- function(input, output, session, rv) {
    observe({
    lapply(rv$table_names, function(tbl_name) {

        raw_heatmap_plot <- ExtendedTask$new(function(raw_input){
          future_promise({
            raw_data <- raw_input$raw_data
            timepoint_cols <- raw_input$timepoint_cols
            datatype <- raw_input$datatype

            clean_data <- as.matrix(raw_data[timepoint_cols])
            if (datatype == 'proteomics') {
              rownames(clean_data) <- raw_data$Gene_Name
            } else if (datatype == 'rnaseq') {
              rownames(clean_data) <- raw_data$Gene_Name
            } else if (datatype == 'phosphoproteomics') {
              rownames(clean_data) <- raw_data$pepG
            }
            raw_ht <- Heatmap(clean_data, show_row_dend = F, show_row_names = F)

            list(raw_ht = raw_ht)
          })
        })

        observeEvent(input[[paste0("compute_raw_ht_", tbl_name)]], {
  
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
          
        output[[paste0("raw_ht_", tbl_name)]] <- renderPlot({
          result <- raw_heatmap_plot$result()
          raw_ht <- result$raw_ht
          raw_ht
        })
  
        
      
    })
  })

}


#### THIS FUNCTION CONTAINS BOTH THE SERVER FOR THE DEP HEATMAP AND THE VOLCANO (because both use the dep_output)

dep_heatmap_server <- function(input, output, session, rv, cache) {
  observe({
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

      req(input[[paste0("volcano_pcutoff_", tbl_name)]],  input[[paste0("volcano_fccutoff_", tbl_name)]])
      pcut <- input[[paste0("volcano_pcutoff_", tbl_name)]]
      lfcut <- input[[paste0("volcano_fccutoff_", tbl_name)]]

      dep_output <- add_rejections(dep_output, alpha = 10^-pcut, lfc = lfcut)
      rd_names <- colnames(rowData(dep_output))
      sig_cols <- grep("__significant$", rd_names, value = TRUE)
      valid_contrasts <- sub("__significant$", "_", sig_cols)  
      #Isolate function to avoid retriggering the observe block
      isolate({ 
        rv$dep_output[[tbl_name]] <- dep_output
        rv$contrasts[[tbl_name]] <- valid_contrasts
      })  

#### TESTING NON BLOCKING

     volcano_plot <- ExtendedTask$new(function(deps) {
        future_promise({

          dep_output <- deps$dep_output
          contrast <- deps$contrast
          datatype <- deps$datatype
          pcut <- deps$pcut
          fccut <- deps$fccut
          highlight <- deps$highlight

          lfc_col <- paste0(contrast, "_diff")
          pval_col <- paste0(contrast, "_p.adj")
          res <- rowData(dep_output)

          # Prepare base values
          df_table <- NULL
          labels <- NULL
          text_col <- NULL

          if (datatype == 'proteomics') {
            gene_names <- str_to_title(res$Gene_Name)
            log2FC <- res[[lfc_col]]
            pval <- res[[pval_col]]
            df <- data.frame(gene_names, pval, log2FC)
            sig_col <- paste0(contrast, "_significant")
            if (all(is.na(sig_col)) || all (sig_col == FALSE)) {
              showNotification("No significant genes for this contrast.", type = "error")
              return()
            } else {
              df_table <- get_results(get_signicant(dep_output, contrast))
            }
            text_col <- gene_names

            volc <- EnhancedVolcano(
              df,
              lab = gene_names,
              selectLab = c("a"),
              x = 'log2FC',
              y = 'pval',
              title = "",
              pCutoff = pcut,
              FCcutoff = fccut,
              pointSize = ifelse(text_col %in% highlight, 3, 1),
              legendPosition = 'none'
            )
            volc <- volc + aes(text = gene_names) + labs(color = "Legend")

          } else if (datatype == 'phosphoproteomics') {
            gene_names <- str_to_title(res$Gene_Name)
            peptide <- res$pepG
            log2FC <- res[[lfc_col]]
            pval <- res[[pval_col]]

            df <- data.frame(peptide, pval, log2FC)

            sig_se <- get_signicant(dep_output, contrast)
            df_table <- get_results(sig_se)
            res_all <- rowData(dep_output)
            res_sig <- rowData(sig_se)
            sig_genes <- str_to_title(res_sig$Gene_Name)
            df_table$Gene_Name <- str_to_title(res_sig$Gene_Name)
            df_table <- df_table[, c("Gene_Name", colnames(df_table)[colnames(df_table) != "Gene_Name"])]

            text_col <- peptide

            volc <- EnhancedVolcano(
              df,
              lab = peptide,
              selectLab = c("a"),
              x = 'log2FC',
              y = 'pval',
              title = "",
              pCutoff = pcut,
              FCcutoff = fccut,
              pointSize = ifelse(text_col %in% highlight, 3, 1),
              legendPosition = 'none'
            )
            volc <- volc + aes(text = peptide) + labs(color = "Legend")

          } else if (datatype == 'rnaseq') {
            gene_ID <- rownames(res)
            log2FC <- res[[lfc_col]]
            pval <- res[[pval_col]]
            df <- data.frame(gene_ID, pval, log2FC)

            df_table <- get_results(get_signicant(dep_output, contrast))
            text_col <- gene_ID

            volc <- EnhancedVolcano(
              df,
              lab = gene_ID,
              selectLab = c("a"),
              x = 'log2FC',
              y = 'pval',
              title = "",
              pCutoff = pcut,
              FCcutoff = fccut,
              legendPosition = 'none'
            )
            volc <- volc + labs(color = "Legend")
          }


          list(volcano = volc, table = df_table, df = df)
        })
      })

      observeEvent(input[[paste0("compute_volcano_", tbl_name)]], {
        contrast <- isolate(input[[paste0("comparison_volcano_", tbl_name)]])
        pcut <- isolate(10^(-input[[paste0("volcano_pcutoff_", tbl_name)]]))
        fccut <- isolate(input[[paste0("volcano_fccutoff_", tbl_name)]])
        highlight <- isolate(input[[paste0("volcano_search_", tbl_name)]]) %>%
          str_split(",", simplify = TRUE) %>%
          str_trim() %>%
          str_to_title()
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


      output[[paste0("volcano_", tbl_name)]] <- renderPlotly({
        result <- volcano_plot$result()
        req(result)
        df <- result$df
        ggplotly(result$volcano + aes(x = log2FC, y = -log10(pval)), tooltip = "text")
      })

      output[[paste0("volcano_sig_table_", tbl_name)]] <- renderDT({
        result <- volcano_plot$result()
        req(result)
        datatable(result$table, options = list(scrollX = TRUE, pageLength = 10))
      })



#### END OF TESTING NON BLOCKING

      plot_dep_heatmap <- ExtendedTask$new(function(ht_inps) {
          future_promise({
            clustering_enabled <- ht_inps$clustering_enabled
            num_clusters <- ht_inps$num_clusters
            dep_output <- ht_inps$dep_output
            if (clustering_enabled) {
            # Apply clustering based on num_clusters
              dep_output_plot <- plot_heatmap(dep_output, kmeans = TRUE, k = num_clusters)
            } else {
              dep_output_plot <- plot_heatmap(dep_output)
            }
            dep_pg_sig <- get_signicant(dep_output)
            expr <- assay(dep_pg_sig)
            gene_info <- as.data.frame(rowData(dep_pg_sig))
            df <- cbind(gene_info, as.data.frame(expr))
            df <- df[, c(colnames(gene_info), colnames(expr))]
            list(heatmap = dep_output_plot, df = df)
          })
        })

        observeEvent(input[[paste0("recompute_heatmap_", tbl_name)]], {
          # Check clustering option
          clustering_enabled <- isolate(input[[paste0("clustering_", tbl_name)]])
          num_clusters <- isolate(input[[paste0("num_clusters_", tbl_name)]])
          dep_output <- isolate(rv$dep_output[[tbl_name]])

          ht_inps <- list(
            clustering_enabled = clustering_enabled,
            num_clusters = num_clusters,
            dep_output = dep_output
          )

          plot_dep_heatmap$invoke(ht_inps)
        })
  
 
        output[[paste0("ht_", tbl_name)]] <- renderPlot({
          result <- plot_dep_heatmap$result()
          dep_output_plot <- result$heatmap
          dep_output_plot
        })

        output[[paste0("ht_sig", tbl_name)]] <- renderDT({
          result <- plot_dep_heatmap$result()
          df <- result$df
          datatable(df, options = list(scrollX = T, pageLength = 10))
        })
      
    })
  })
}


