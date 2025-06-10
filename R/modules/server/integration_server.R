source("R/global.R")
source("R/modules/server/int_timeline_server.R")
source("R/modules/helpers/integration_helper.R")

integration_ui <- function(input, output, session, rv) {
  combined_data <- reactiveVal(NULL)

  output$integration_ui <- renderUI({
    req(length(rv$tables) > 0)
    tagList(
        tabBox(title = "Integration", id = "integration_tabs", selected = "Raw Integration", width = 12, 
            tabPanel("Raw Integration",
                fluidRow(
                    box(title = "Integration Settings", width = 12, status = "success", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                    fluidRow(
                        box(
                        title = "Integrate", width = 4, solidHeader = TRUE, status = "primary",
                        pickerInput(inputId = "integration", label = "Integrate", choices = rv$table_names, multiple = TRUE, options = pickerOptions(container = "body"), width = "100%"),
                        uiOutput("integration_col_selector"),
                        div(style = "text-align: center;",
                            actionBttn("integrate_data", span("Integrate", style = "color: black;"), icon = span(icon("arrow-right-to-bracket"), style = "color: black;"), style = "jelly", size = "sm", class = "btn-primary")
                        )
                        ),
                        uiOutput("preview_box")
                    )
                ), 
                box(
                    title = "Integrated Raw Table", width = 12, solidHeader = TRUE, status = "info", collapsible = TRUE, collapsed = FALSE, 
                    DTOutput("integration_combined_table")
                ),
                box(
                    title = "Integrated Timeline Plot", width = 12, solidHeader = TRUE, status = "info",
                    selectizeInput(
                        inputId = "search_gene_integration",
                        label = "Search your gene of interest:",
                        choices = NULL,
                        multiple = TRUE
                    ),
                    selectInput(
                        inputId = "scale_integration",
                        label = "Select scale:",
                        choices = c("Continous", "Log-scale"),
                        selected = "Continous"
                    ),
                    plotOutput("integration_timeline_plot")
                )
            )
            ),
            tabPanel("Processed Integration",
                fluidRow(
                    box(title = "Processed Integration Settings", width = 12, status = "success", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                    fluidRow(
                        box(
                            title = "Process Integrate", width = 4, solidHeader = TRUE, status = "primary",
                            pickerInput(
                                inputId = "processed_integration", label = "Integrate", 
                                choices = rv$table_names, multiple = TRUE, 
                                options = pickerOptions(container = "body"), width = "100%"
                            ),
                            p("Select for each table upon which contrast do the filtering (they should match)"),
                            uiOutput("comparison_col_selector_pi"),
                            p("----------"),
                            numericInput("heatmap_k", "Number of clusters (k):", value = 3, min = 2, max = 10),
                            numericInput("pval_thresh_pi", "P-value threshold:", value = 0.05, min = 0, step = 0.01),
                            numericInput("lfc_thresh_pi", "LFC threshold:", value = 1, min = 0, step = 0.1),

                            div(style = "text-align: center;",
                                actionBttn("process_integrate_data", span("Integrate", style = "color: black;"), 
                                            icon = span(icon("arrow-right-to-bracket"), style = "color: black;"), 
                                            style = "jelly", size = "sm", class = "btn-primary")
                            )
                            ),
                            uiOutput("preview_box_integrate")
                    )),
                    box(title = "Tables", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, 
                        fluidRow(
                            uiOutput("processed_integrated_tables")
                        )
                    ),
                    box(title = "Integrated Heatmaps", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                      uiOutput("integrated_heatmaps")
                    ),
                    box(title = "LFC Scatterplot", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                      uiOutput("lfc_scatter_ui")
                    )

                ))
       )
    )
  })

  # Preview box with column matching and unified table
  output$preview_box <- renderUI({
    req(length(input$integration) > 0)
    tagList(
      box(
        title = "Integration Preview", width = 8, solidHeader = TRUE, status = "primary",
        uiOutput("integration_column_matching")
      )
    )
  })



  output$comparison_col_selector_pi <- renderUI({
    req(input$processed_integration)
    lapply(input$processed_integration, function(tbl) {
      available_comparisons <- rv$contrasts[[tbl]]
      selectInput(
        inputId = paste0("pi_comparison_selected_", tbl),
        label = paste0("Select comparison for ", tbl),
        choices = available_comparisons,
        selected = NULL,
        width = "100%", 
        multiple = FALSE
      )
    })
  })



    output$preview_box_integrate <- renderUI({
      req(rv$integration_preview_dims)

      dims_ui <- lapply(names(rv$integration_preview_dims), function(tbl) {
        dims <- rv$integration_preview_dims[[tbl]]
        tagList(
          tags$li(tags$b(tbl)),
          tags$ul(
            tags$li(paste("Original:", paste(dims$original, collapse = " x "))),
            tags$li(paste("After filtering:", paste(dims$filtered, collapse = " x "))),
            tags$li(paste("After intersection:", paste(dims$intersected, collapse = " x ")))
          )
        )
      })

      tagList(
        box(
          title = "Processed Integration Preview", width = 8, solidHeader = TRUE, status = "primary",
          tags$ul(
            tags$li(tags$b("Table dimensions at each step:")),
            dims_ui,
            tags$li(tags$b("Selected comparisons:")),
            tags$ul(
              lapply(input$processed_integration, function(tbl) {
                comp <- input[[paste0("pi_comparison_selected_", tbl)]]
                tags$li(paste(tbl, ":", comp))
              })
            )
          )
        )
      )
    })


    output$processed_integrated_tables <- renderUI({
        req(rv$intersected_tables_processed)

        tables <- rv$intersected_tables_processed

        # Generate a UI list of DTOutput placeholders
        tagList(
            lapply(names(tables), function(tbl) {
            box(
                title = paste("Processed Table:", tbl),
                width = 12,
                solidHeader = TRUE,
                status = "info",
                DTOutput(outputId = paste0("processed_tbl_", tbl))
            )
            })
        )
    })

    observe({
        req(rv$intersected_tables_processed)
        lapply(names(rv$intersected_tables_processed), function(tbl) {
            local({
            table_name <- tbl
            output_id <- paste0("processed_tbl_", table_name)
            data <- as.data.frame(rv$intersected_tables_processed[[table_name]])
            output[[output_id]] <- renderDT({
                datatable(data, options = list(scrollX = TRUE, pageLength = 5))
            })
            })
        })
    })




  # Column selector for each selected table
  output$integration_col_selector <- renderUI({
    req(length(input$integration) > 0)
    lapply(input$integration, function(int_tbl) {
      cols <- rv$time_cols[[int_tbl]]
      pickerInput(
        inputId = paste0("cols_selected_int", int_tbl),
        label = paste0("Select columns from ", int_tbl, " to load:"),
        choices = c(cols),
        selected = NULL,
        multiple = TRUE,
        options = pickerOptions(
          actionsBox = TRUE,
          liveSearch = TRUE,
          noneSelectedText = "Select columns",
          deselectAllText = "None",
          selectAllText = "Select All",
          dropupAuto = FALSE
        ),
        width = "100%"
      )
    })
  })




  # Show how columns will be matched
  output$integration_column_matching <- renderUI({
    req(length(input$integration) > 1)

    selected_columns <- lapply(input$integration, function(tbl) {
      input[[paste0("cols_selected_int", tbl)]]
    })

    if (length(unique(sapply(selected_columns, length))) > 1) {
      return(p("Column numbers differ. Matching not possible.", style = "color: red;"))
    }

    n <- length(selected_columns[[1]])
    if (n == 0) return(NULL)

    matches <- lapply(seq_len(n), function(i) {
      paste(
        sapply(seq_along(selected_columns), function(j) {
          selected_columns[[j]][i]
        }),
        collapse = " ⇄ "
      )
    })

    tagList(
      tags$ul(
        lapply(matches, function(line){
          tags$li(line)
        })
      )
    )
  })

  # Combined joined table based on selected columns and shared identifiers
  observeEvent(input$integrate_data, {
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


####### OBSERVE EVENT PROCESSED DATA ######

    observeEvent(input$process_integrate_data, {
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
        
        # Extract the necessary columns from rowData
        res <- as.data.frame(rowData(dep))

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
      
        original_dim <- dim(assay(dep))
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

      # Subset assays by intersected gene names
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
      }

      # Save results
      rv$intersected_tables_processed <- intersected_list
      rv$integration_preview_dims <- lapply(selected_tables, function(tbl) {
        list(
          original = dim(assay(rv$dep_output[[tbl]])),
          filtered = dim_info[[tbl]]$filtered,
          intersected = dim_info[[tbl]]$intersected
        )
      })
      names(rv$integration_preview_dims) <- selected_tables
    })


  output$integrated_heatmaps <- renderUI({
    req(rv$intersected_tables_processed, input$heatmap_k)
    tables <- rv$intersected_tables_processed
    tagList(
        lapply(names(rv$intersected_tables_processed), function(tbl) {
        plotOutput(outputId = paste0("heatmap_", tbl), width = "80%")
        }),
        br(),
        lapply(names(tables), function(tbl) {
          box(
            title = paste("Clusters Table:", tbl),
            width = 12,
            solidHeader = TRUE,
            status = "info",
            collapsible = TRUE, 
            collapsed = FALSE, 
            DTOutput(outputId = paste0("cluster_table_", tbl))
          )
        })
    )   
  })

  observe({
    req(rv$intersected_tables_processed, input$heatmap_k)

    all_tables <- rv$intersected_tables_processed
    k <- input$heatmap_k

    first_tbl <- all_tables[[1]]

    mat_scaled <- t(scale(t(first_tbl)))  # Scale across rows

    good_rows <- apply(mat_scaled, 1, function(x) all(!is.na(x) & !is.nan(x)))
    mat_scaled <- mat_scaled[good_rows, ]

    # Perform k-means clustering on scaled data
    if (nrow(mat_scaled) < k) {
      showNotification(paste("Cannot create", k, 
                            "clusters from", nrow(mat_scaled), "genes."),
                      type = "error")
      return(NULL)
    }
    km <- kmeans(mat_scaled, centers = k)

    # Get cluster labels and order of genes
    cluster_labels <- km$cluster
    ordered_genes <- rownames(mat_scaled)[order(cluster_labels)]

    # Save info for use across plots
    rv$heatmap_clusters <- list(
        order = ordered_genes,
        cluster = cluster_labels
    )

    gene_map <- data.frame(unique_id = rownames(first_tbl), Gene_Name = sapply(strsplit(rownames(first_tbl), "_"), `[`, 1))
    rv$gene_to_cluster <- setNames(rv$heatmap_clusters$cluster, gene_map$Gene_Name)

    # Render heatmaps
    lapply(names(all_tables), function(tbl) {
        local({
        tbl_name <- tbl
        mat <- all_tables[[tbl_name]]
        mat_scaled_tbl <- t(scale(t(mat)))

        gene_names <- sapply(strsplit(rownames(mat_scaled_tbl), "_"), `[`, 1)
        cluster_vec <- factor(rv$gene_to_cluster[gene_names], levels = 1:k)

        valid <- !is.na(cluster_vec)
        mat_ordered <- mat_scaled_tbl[valid, , drop = FALSE]
        cluster_vec <- cluster_vec[valid]

        ordered <- order(cluster_vec)
        mat_ordered <- mat_ordered[ordered, , drop = FALSE]
        cluster_vec <- cluster_vec[ordered]

        output[[paste0("heatmap_", tbl_name)]] <- renderPlot({
          mat <- all_tables[[tbl_name]]
          mat_scaled_tbl <- t(scale(t(mat)))

          # Build cluster average profiles
          cluster_ids <- split(seq_len(nrow(mat_ordered)), cluster_vec)
          line_profiles <- t(vapply(cluster_ids, function(idxs) {
            colMeans(mat_ordered[idxs, , drop = FALSE], na.rm = TRUE)
          }, FUN.VALUE = numeric(ncol(mat_ordered))))

          # Safeguard: avoid empty trend lines
          if (nrow(line_profiles) == 0 || ncol(line_profiles) == 0) {
            showNotification("Could not compute trend lines — data missing or clustering failed", type = "error")
            return(NULL)
          }

          # Normalize trend lines for plotting
          line_profiles_norm <- t(apply(line_profiles, 1, function(x) {
            rng <- range(x, na.rm = TRUE)
            if (diff(rng) == 0) rep(0.5, length(x)) else (x - rng[1]) / diff(rng)
          }))

          # Annotation: line plots next to each cluster
          trend_anno <- rowAnnotation(trend = anno_link(
            align_to = cluster_vec,
            which = "row",
            panel_fun = function(index, nm) {
              grid.rect()
              grid.lines(
                x = seq_len(ncol(line_profiles_norm)) / ncol(line_profiles_norm),
                y = line_profiles_norm[as.integer(nm), ],
                gp = gpar(col = "#2b8cbe", lwd = 1)
              )
            },
            side = "right",
            size = unit(3, "cm"),
            width = unit(5, "cm")
          ))
          # Final heatmap with trend annotation
          ComplexHeatmap::Heatmap(
            mat_ordered,
            name = tbl_name,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = FALSE,
            show_column_dend = TRUE,
            row_split = cluster_vec,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 8),
            heatmap_legend_param = list(title = "Expression"),
            right_annotation = trend_anno
          )
        })


      })
    })

    lapply(names(all_tables), function(tbl) {
      local({
        tbl_name <- tbl
        data <- all_tables[[tbl_name]]

        gene_names <- sapply(strsplit(rownames(data), "_"), `[`, 1)
        clusters <- rv$gene_to_cluster[gene_names]
        clusters[is.na(clusters)] <- NA

        df <- data.frame(
          Unique_ID = rownames(data),
          Gene_Name = gene_names,
          Cluster = clusters,
          stringsAsFactors = FALSE
        )

        output[[paste0("cluster_table_", tbl_name)]] <- DT::renderDT({
          DT::datatable(df, options = list(scrollX = TRUE, pageLength = 5))
        })
      })
    })


    output$lfc_scatter_ui <- renderUI({
      req(rv$intersected_tables_processed)
      selected_tables <- names(rv$intersected_tables_processed)
      req(length(selected_tables) == 2)

      if (any(rv$datatype[selected_tables] == "phosphoproteomics")) {
        # Return a nicely styled message div
        div(
          style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
          "Scatter plot is not available when phosphoproteomics data is included because of non-1:1 row mapping."
        )
      } else {
        # Return the plotly output UI placeholder
        plotlyOutput("lfc_scatter_plot", height = "400px", width = "100%")
      }
    })



    output$lfc_scatter_plot <- renderPlotly({
      req(rv$intersected_tables_processed)
      selected_tables <- names(rv$intersected_tables_processed)
      req(length(selected_tables) == 2)


      scatter_data <- data.frame(Gene_Name = rownames(rv$intersected_tables_processed[[1]]))

      for (tbl in selected_tables) {
        contrast <- input[[paste0("pi_comparison_selected_", tbl)]]
        dep <- rv$dep_output[[tbl]]
        lfc_col <- paste0(contrast, "_diff")

        res <- as.data.frame(rowData(dep))

        if (rv$datatype[[tbl]] == "rnaseq"){
          mapper <- rv$tables[[tbl]]
          mapper <- mapper[, c("Gene_ID", "Gene_Name")]
          res$Gene_ID <- rownames(res)
          res <- merge(res, mapper, by = "Gene_ID", all.x = TRUE)
        }

        df <- res[, c("Gene_Name", lfc_col)]
        colnames(df) <- c("Gene_Name", paste0("LFC: ", tbl))
        scatter_data <- merge(scatter_data, df, by = "Gene_Name", all.x = TRUE)
      }

      colnames_data <- colnames(scatter_data)

      # Plot using ggplot and plotly
      p <- ggplot(scatter_data, aes(x = .data[[colnames_data[2]]], y = .data[[colnames_data[3]]], text = Gene_Name)) +
        geom_point(alpha = 0.7, color = "#2b8cbe") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        labs(
          title = "LFC Scatter Plot of Intersected Genes"
        ) +
        theme_minimal()

      ggplotly(p, tooltip = "text")
    })


  })




  # Render integrated table
    output$integration_combined_table <- renderDT({
    if (is.null(combined_data())) {
        return(
        datatable(
            data.frame(Message = "Click Integrate button to see the Integrated Table"),
            options = list(
            dom = 't',
            ordering = FALSE,
            paging = FALSE,
            searching = FALSE
            ),
            rownames = FALSE
        )
        )
    }

    DT::datatable(combined_data(), options = list(scrollX = TRUE, pageLength = 10))
    })



}
