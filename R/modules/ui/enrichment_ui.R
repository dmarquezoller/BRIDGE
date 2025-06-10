enrichment_settings_ui <- function(tbl_name, rv) {
    box(title = "Settings", width = 5, solidHeader = TRUE, status = "info", 
        virtualSelectInput(paste0("comparison_db_", tbl_name), "Select where to perform enrichment analysis:", 
            choices = c("GO", "KEGG", "REACTOME"),
            selected = "GO"),
        virtualSelectInput(paste0("contrasts_enrichment_", tbl_name), "Select Contrast:", 
            choices = rv$contrasts[[tbl_name]], 
            selected = rv$contrasts[[tbl_name]][1]
        ),
        actionBttn(paste0("compute_enrichment_", tbl_name), span("Compute Enrichment Plot", style = "color: white;"), style = "simple", color = "primary", size = "sm")

                           
    )
}

enrichment_plots_ui <- function(tbl_name) {
    box(title = "Enrichment plots", width = 7, solidHeader = T, status = "info",
        h3(),
        uiOutput(paste0("enrichment_", tbl_name)), 
    )
}