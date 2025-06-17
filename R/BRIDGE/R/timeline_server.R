#source("R/global.R")
#' @export
timeline_server <- function(input, output, session, rv) {
    observe({
    lapply(rv$table_names, function(tbl_name) {
  
    if (rv$datatype[[tbl_name]] == "phosphoproteomics") {
      gene <- str_to_lower(trimws(input[[paste0("search_gene_", tbl_name)]]))
      raw_data <- rv$tables[[tbl_name]]
      req(input[[paste0("scale_", tbl_name)]])
      scale_input <- input[[paste0("scale_", tbl_name)]]
      if ( scale_input == "Total Intensity") {
        processed_data <- total_intensity_normalization(raw_data)
      } else if (scale_input == "Continous") {
        processed_data <- raw_data
      } else if (scale_input == "Log-scale") {
        processed_data <- log2_transform(raw_data)
      } else if (scale_input == "Median Normalization") {
        processed_data <- median_normalization(raw_data)
      } else if (scale_input == "FPKM") {
        annotation <- rv$annotation[[tbl_name]]
        processed_data <- fpkm_normalization(raw_data, annotation)
      } else if (scale_input == "TMM") {
        processed_data <- tmm_normalization(raw_data)
      } else if (scale_input == "CPM") {
        processed_data <- cpm_normalization(raw_data)
      } else if (scale_input == "TPM") {
        annotation <- rv$annotation[[tbl_name]]
        fpkm_data <- fpkm_normalization(raw_data, annotation)
        processed_data <- tpm_normalization(fpkm_data)
      }
      data <- subset(processed_data, str_to_lower(Gene_Name) %in% gene)
      data$Gene_pepG <- paste(str_to_title(data$Gene_Name), data$pepG, sep = "_")
      
      timepoints <- c("Gene_pepG", rv$time_cols[[tbl_name]])
      unique_timepoints <- rv$time_cols[[tbl_name]] %>%
        gsub("_[0-9]+$", "", .) %>%  # Remove trailing _1, _2, etc.
        gsub("[_.-]+$", "", .) %>%   # Clean up any extra underscores/dots
        unique() 
      data <- data[ ,timepoints]
      
      validate(need(nrow(data) > 0, "No data found for the entered gene(s)"))
      data_long <- data %>%
        pivot_longer(cols=-Gene_pepG, names_to = "Stage", values_to = "Expression") %>%
        mutate(StageGroup = gsub("[0-9]+$", "", Stage) %>% gsub("[_.-]+$", "", .))
      data_long$StageGroup <- factor(data_long$StageGroup, 
                                     levels = unique_timepoints)
      data_avg <- data_long %>%
        group_by(StageGroup, Gene_pepG) %>%
        summarize(MeanExpression = mean(Expression, na.rm=T), .groups = "drop")
      
      output[[paste0("time_plot_dt_", tbl_name)]] <- DT::renderDT({
        datatable(data)
      })
      
      
      output[[paste0("time_plot_", tbl_name)]] <- renderPlot({
        
        plot <- ggplot() +
          # Plot individual replicates as dots (aligned vertically)
          geom_point(data = data_long, 
                     aes(x = StageGroup, y = Expression, color = as.factor(str_to_title(Gene_pepG))), 
                     size = 3, alpha = 0.6) +  
          
          # Plot the mean expression line
          geom_line(data = data_avg, 
                    aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Gene_pepG)), group = str_to_title(Gene_pepG)), 
                    linewidth = 1.2) +
          
          # Add mean expression points (triangles)
          geom_point(data = data_avg, 
                     aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Gene_pepG))), 
                     size = 4, shape = 17) +
          # Labels and theme
          labs(
            x = "Stage",
            y = "log Expression",
            title = paste("Peptide Expression Time Line", str_to_title(input[[paste0("search_gene_", tbl_name)]])),
            color = "Proteins"
          ) +
          theme_minimal() +
          theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust= 1)) 
        if (input[[paste0("scale_", tbl_name)]] == "Log-scale"){
          plot <- plot + scale_y_log10()
          
        }
        return(plot)

        
      })
    }else{
        gene <- str_to_lower(trimws(input[[paste0("search_gene_", tbl_name)]]))
        raw_data <- rv$tables[[tbl_name]]
        req(input[[paste0("scale_", tbl_name)]])
        scale_input <- input[[paste0("scale_", tbl_name)]]
        if ( scale_input == "Total Intensity") {
          processed_data <- total_intensity_normalization(raw_data)
        } else if (scale_input == "Continous") {
          processed_data <- raw_data
        } else if (scale_input == "Log-scale") {
          processed_data <- log2_transform(raw_data)
        } else if (scale_input == "Median Normalization") {
          processed_data <- median_normalization(raw_data)
        } else if (scale_input == "FPKM") {
          annotation <- rv$annotation[[tbl_name]]
          processed_data <- fpkm_normalization(raw_data, annotation)
        } else if (scale_input == "TMM") {
          processed_data <- tmm_normalization(raw_data)
        } else if (scale_input == "CPM") {
          processed_data <- cpm_normalization(raw_data)
        } else if (scale_input == "TPM") {
          annotation <- rv$annotation[[tbl_name]]
          fpkm_data <- fpkm_normalization(raw_data, annotation)
          processed_data <- tpm_normalization(fpkm_data)
        }
        data <- subset(processed_data, str_to_lower(Gene_Name) %in% gene)
        timepoints <- c("Gene_Name", rv$time_cols[[tbl_name]])
        unique_timepoints <- rv$time_cols[[tbl_name]] %>%
          gsub("_[0-9]+$", "", .) %>%  # Remove trailing _1, _2, etc.
          gsub("[_.-]+$", "", .) %>%   # Clean up any extra underscores/dots
          unique() 
        data <- data[ ,timepoints]
        validate(need(nrow(data) > 0, "No data found for the entered gene(s)"))

        
        data_long <- data %>%
          pivot_longer(cols = -Gene_Name, names_to = "Stage", values_to = "Expression") %>%
          mutate(StageGroup = gsub("[0-9]+$", "", Stage) %>% gsub("[_.-]+$", "", .))
        data_long$StageGroup <- factor(data_long$StageGroup, 
                                       levels = unique_timepoints)
        data_avg <- data_long %>%
          group_by(StageGroup, Gene_Name) %>%
          summarize(MeanExpression = mean(Expression, na.rm=T), .groups = "drop")
        
        output[[paste0("time_plot_dt_", tbl_name)]] <- DT::renderDT({
          datatable(data)
        })
        
        output[[paste0("time_plot_",tbl_name)]] <- renderPlot({
  
            plot <- ggplot() +
              # Plot individual replicates as dots (aligned vertically)
              geom_point(data = data_long, 
                         aes(x = StageGroup, y = Expression, color = as.factor(str_to_title(Gene_Name))), 
                         size = 3, alpha = 0.6) +  
              
              # Plot the mean expression line
              geom_line(data = data_avg, 
                        aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Gene_Name)), group = str_to_title(Gene_Name)), 
                        linewidth = 1.2) +
              
              # Add mean expression points (triangles)
              geom_point(data = data_avg, 
                         aes(x = StageGroup, y = MeanExpression, color = as.factor(str_to_title(Gene_Name))), 
                         size = 4, shape = 17) +
              # Labels and theme
              labs(
                x = "Stage",
                y = "log Expression",
                title = paste("Protein Expression Time Line", str_to_title(input[[paste0("search_gene_", tbl_name)]])),
                color = "Proteins"
              ) +
              theme_minimal() +
              theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust= 1)) 

              
            return(plot)
        })
    }
    })
  })
}