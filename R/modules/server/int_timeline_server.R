source("R/global.R")

int_timeline_server <- function(input, output, session, data_combined, reference_time_names) {

  processed_data <- reactive({
    req(input$scale_integration)
    scale_input <- input$scale_integration
    raw_data <- data_combined
    if (scale_input == "Continous") {
      raw_data
    } else if (scale_input == "Log-scale") {
      log2_transform(raw_data)
    } 
  })

  selected_data <- reactive({
    req(input$search_gene_integration)
    gene <- str_to_lower(trimws(input$search_gene_integration))
    data <- processed_data()
    data <- subset(data, str_to_lower(Gene_Name) %in% gene)
    validate(need(nrow(data) > 0, "No data found for the entered gene(s)"))
    data
  })

  data_long <- reactive({
    data <- selected_data()

    # Clean stage names
    unique_timepoints <- reference_time_names %>%
      gsub("_[0-9]+$", "", .) %>%
      gsub("[_.-]+$", "", .) %>%
      unique()

    # Reshape
    data_long <- data %>%
      pivot_longer(cols = all_of(reference_time_names), names_to = "Stage", values_to = "Expression") %>%
      mutate(
        StageGroup = gsub("[0-9]+$", "", Stage) %>% gsub("[_.-]+$", "", .),
        StageGroup = factor(StageGroup, levels = unique_timepoints)
      )

    data_long
  })
  data_avg <- reactive({
    data_long() %>%
      group_by(source, StageGroup, unique_id) %>%
      summarize(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
  })

  output$integration_timeline_plot <- renderPlot({
    long <- data_long()
    avg <- data_avg()

    ggplot() +
      geom_point(data = long, aes(x = StageGroup, y = Expression, color = unique_id), size = 3, alpha = 0.5) +
      geom_line(data = avg, aes(x = StageGroup, y = MeanExpression, color = unique_id, group = unique_id), linewidth = 1.2) +
      geom_point(data = avg, aes(x = StageGroup, y = MeanExpression, color = unique_id), size = 4, shape = 17) +
      facet_wrap(~ source, ncol = 2, scales = "free_y") +
      labs(x = "Stage", y = "Expression", color = "Gene") +
      theme_minimal() +
      theme(
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold")
      )
  })

}
