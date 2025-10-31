#' @export
EnrichmentServer <- function(id, rv, tbl_name) {
    moduleServer(id, function(input, output, session) {
        # Populate contrast choices when available
        strip_sig <- function(dep_obj) {
            if (methods::is(dep_obj, "DEGdata")) {
                tr <- dep_obj@test_result
                tr <- tr[, !grepl("_significant$|^significant$", colnames(tr))]
                dep_obj@test_result <- tr
            } else {
                rd <- SummarizedExperiment::rowData(dep_obj)
                rd <- rd[, !grepl("_significant$|^significant$", colnames(rd))]
                SummarizedExperiment::rowData(dep_obj) <- rd
            }
            dep_obj
        }

        check_sig <- function(dep_obj) {
            if (methods::is(dep_obj, "DEGdata")) {
                tr <- as.data.frame(dep_obj@test_result)
                check <- nrow(tr %>% dplyr::filter(significant == TRUE))
            } else {
                rd <- as.data.frame(SummarizedExperiment::rowData(dep_obj))
                check <- nrow(rd %>% dplyr::filter(significant == TRUE))
            }
            # message("Significant hits: ", check)
            check
        }

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
                DEP2::check_organismDB_depends(species, install = TRUE)                
                # Fallback for unsupported organisms
                Gene_Names <- stringr::str_to_lower(gsub("_.*$", "", rownames(dep_output)))
                #rownames(dep_output) <- stringr::str_to_upper(rownames(dep_output))
                rownames(dep_output) <- Gene_Names
                #message("NAMES: ", paste(rownames(dep_output)[1:5], collapse = ", "))

                if (!(species %in% c("Zebrafish", "Human", "Mouse"))) {                 
                    species <- "Human"
                }

                dep_output <- strip_sig(dep_output) # Remove old sig cols if present
                # Sig filtering used in your original code before ORA
                dep_pg <- DEP2::add_rejections(dep_output, alpha = input$enrichment_pcutoff, lfc = input$enrichment_fccutoff)                
                sig_count <- check_sig(dep_pg)
                #message("Significant hits for enrichment: ", sig_count)
                if (sig_count < 10) {                
                    output$enrichment <- renderUI({
                        div(
                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                            paste0("Not enough significant hits (n < 10) for enrichment analysis with contrast '", contrast, "'.")
                        )
                    })
                    return(invisible())
                } else{
                    message("Performing ORA enrichment analysis with contrast '", contrast, "' species '", species, "' and type '", enrichment_type, "' on '", sig_count, "' hits like '",paste(rownames(dep_pg)[1:5], collapse = ',') ,"' ...")
                    
                    # Try/catch for DEP2::test_ORA
                    res_ora <- tryCatch({
                        withCallingHandlers(
                            DEP2::test_ORA(
                                dep_pg,
                                type = enrichment_type,
                                species = species,
                                contrasts = contrast
                            ),
                            warning = function(w) {
                                if (grepl("no applicable method for ", conditionMessage(w))) {
                                    output$enrichment <- renderUI({
                                        div(
                                            style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                            "No genes matched with annotation in the selected organism database! Please check your cutoffs or try another contrast."
                                        )
                                    })
                                    invokeRestart("muffleWarning")
                                }
                            }
                        )
                    }, error = function(e) {
                        if (grepl("no applicable method for ", conditionMessage(e))) {
                            output$enrichment <- renderUI({
                                div(
                                    style = "padding: 20px; color: #d9534f; font-weight: bold; text-align: center;",
                                    "No genes matched with annotation in the selected organism database! Please check your cutoffs or try another contrast."
                                )
                            })
                            return(NULL)
                        } else {
                            stop(e)
                        }
                    })

                    if (!is.null(res_ora)) {
                        output$enrichment <- renderUI({
                            plotOutput(session$ns("enrichment_plot"), height = "520px")
                        })
                        output$enrichment_plot <- renderPlot({
                            enrichplot::dotplot(res_ora)
                        })
                    }
                }
            },
            ignoreInit = TRUE
        )
    })
}
