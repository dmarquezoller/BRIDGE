source("R/global.R")


enrichment_server <- function(input, output, session, rv) {
    observe(
        lapply(rv$table_names, function(tbl_name) {
            observeEvent(input[[paste0("compute_enrichment_", tbl_name)]], {
                if (rv$datatype[[tbl_name]] == "phosphoproteomics") {
                    output[[paste0("enrichment_", tbl_name)]] <- renderUI({
                        tags$div(
                        style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                        "No enrichment plot available for phosphoproteomics datasets."
                        )
                    })
                    return()  # Skip further processing
                }
                
                species <- str_to_title(rv$species[[tbl_name]])
                enrichment_type <- input[[paste0("comparison_db_", tbl_name)]]
                contrast <- input[[paste0("contrasts_enrichment_", tbl_name)]]
                dep_output <- rv$dep_output[[tbl_name]]
                dep_pg <- dep_output %>% add_rejections(alpha = 0.05,lfc = 2)
                check_enrichment_depends()
                check_organismDB_depends(species) # organism annotation for mouse
                print(species)
                print(enrichment_type)
                print(contrast)
                if (!(species %in% c("Zebrafish", "Human", "Mouse")) ) {
                    rownames(dep_pg) <- str_to_upper(rownames(dep_pg))
                    species = "Human"
                }
                res_ora <- DEP2::test_ORA(dep_pg, type = enrichment_type, species = species, contrasts = contrast)

                output[[paste0("enrichment_", tbl_name)]] <- renderUI({
                    plotOutput(paste0("enrichment_plot_", tbl_name))
                })
                output[[paste0("enrichment_plot_", tbl_name)]] <- renderPlot({
                    enrichplot::dotplot(res_ora)
            })
            })

        })
    )
}