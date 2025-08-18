#' @export
EnrichmentSettingsUI <- function(id, tbl_name, rv) {
    ns <- NS(id)
    shinydashboard::box(title = "Settings", width = 5, solidHeader = TRUE, status = "info", 
        shinyWidgets::virtualSelectInput(ns(paste0("comparison_db_", tbl_name)), "Select where to perform enrichment analysis:", 
            choices = c("GO", "KEGG", "REACTOME"),
            selected = "GO"),
        shinyWidgets::virtualSelectInput(ns(paste0("contrasts_enrichment_", tbl_name)), "Select Contrast:", 
            choices = rv$contrasts[[tbl_name]], 
            selected = rv$contrasts[[tbl_name]][1]
        ),
        shinyWidgets::actionBttn(ns(paste0("compute_enrichment_", tbl_name)), shiny::span("Compute Enrichment Plot", style = "color: white;"), style = "simple", color = "primary", size = "sm")

                           
    )
}

#' @export
EnrichmentPlotsUI <- function(id, tbl_name) {
    ns <- NS(id)
    shinydashboard::box(title = "Enrichment plots", width = 7, solidHeader = T, status = "info",
        h3(),
        shiny::uiOutput(ns(paste0("enrichment_", tbl_name))), 
    )
}