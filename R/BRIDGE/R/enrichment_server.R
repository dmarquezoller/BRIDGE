#' @export
EnrichmentServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        # Populate contrast choices when available
        observe({
            req(rv$contrasts[[tbl_name]])
            updateSelectInput(session, "contrasts_enrichment",
                choices = rv$contrasts[[tbl_name]]
            )
        })

        # Render contrast select UI (so we can safely update it)
        output$contrast_ui <- renderUI({
            selectInput(session$ns("contrasts_enrichment"),
                "Contrast",
                choices = rv$contrasts[[tbl_name]]
            )
        })

        observeEvent(input$compute_enrichment,
            {
                req(rv$dep_output[[tbl_name]], rv$datatype[[tbl_name]], rv$species[[tbl_name]])

                # Phospho: show message and stop
                if (identical(rv$datatype[[tbl_name]], "phosphoproteomics")) {
                    output$enrichment <- renderUI({
                        div(
                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                            "No enrichment plot available for phosphoproteomics datasets."
                        )
                    })
                    return(invisible())
                }

                species <- stringr::str_to_title(rv$species[[tbl_name]])
                enrichment_type <- input$comparison_db
                contrast <- input$contrasts_enrichment
                dep_output <- rv$dep_output[[tbl_name]]

                DEP2::check_enrichment_depends()
                DEP2::check_organismDB_depends(species)

                # Fallback for unsupported organisms
                if (!(species %in% c("Zebrafish", "Human", "Mouse"))) {
                    rownames(dep_output) <- stringr::str_to_upper(rownames(dep_output))
                    species <- "Human"
                }

                # Sig filtering used in your original code before ORA
                dep_pg <- DEP2::add_rejections(dep_output, alpha = 0.05, lfc = 2)

                res_ora <- DEP2::test_ORA(dep_pg,
                    type = enrichment_type,
                    species = species, contrasts = contrast
                )

                output$enrichment <- renderUI({
                    plotOutput(session$ns("enrichment_plot"), height = "520px")
                })
                output$enrichment_plot <- renderPlot({
                    enrichplot::dotplot(res_ora)
                })
            },
            ignoreInit = TRUE
        )
    })
}
