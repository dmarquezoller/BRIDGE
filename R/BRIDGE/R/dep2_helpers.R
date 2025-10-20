#' @export
dep2_proteomics <- function(df, tbl_name, rv) {  
    df <- df %>% mutate(UID = paste0(Gene_Name, "_", Protein_ID))
    unique_pg <- DEP2::make_unique(df, names = "UID", ids = "Protein_ID", delim = ";")
    # unique_pg <- DEP2::make_unique(df, names = "Gene_Name", ids = "Protein_ID", delim = ";")
    ecols <- match(trimws(rv$data_cols[[tbl_name]]), colnames(unique_pg))
    se_pg <- DEP2::make_se_parse(unique_pg, columns = ecols, mode = "delim", sep = "_", remove_prefix = T, log2transform = T)
    # message("SE: ", paste0(colnames(se_pg), collapse = ", "), paste(head(assay(se_pg)), collapse = ", "))
    filter_pg <- DEP2::filter_se(se_pg, thr = 1, fraction = 0.3)
    impute_pg <- DEP2::impute(filter_pg, fun = "MinDet")
    norm_pg <- DEP2::normalize_vsn(impute_pg)
    diff_pg <- DEP2::test_diff(norm_pg, type = "all", fdr.type = "BH")
    sig_pg <- DEP2::add_rejections(diff_pg, alpha = 10^-0.05, lfc = 1)
    return(sig_pg)
}

#' @export
dep2_phosphoproteomics <- function(df, tbl_name, rv) {
    # df <- df[!grepl("p0", df$pepG), ] # Remove the p0   NOT RECOMMENDED!
    df <- df %>% mutate(UID = paste0(Protein_ID, "_", pepG))
    unique_phos <- DEP2::make_unique(df, names = "UID", ids = "UID", delim = ";")
    ecols <- match(trimws(rv$data_cols[[tbl_name]]), colnames(unique_phos))
    # message("Ecols: ", paste0(ecols, collapse = ", "), "\n")
    se_phos <- DEP2::make_se_parse(unique_phos, columns = ecols, mode = "delim", sep = "_", remove_prefix = T, log2transform = T)
    filter_phos <- DEP2::filter_se(se_phos)
    impute_phos <- DEP2::impute(filter_phos, fun = "MinDet")
    norm_phos <- DEP2::normalize_vsn(impute_phos)
    diff_phos <- DEP2::test_diff(norm_phos, type = "all", fdr.type = "BH")
    sig_phos <- DEP2::add_rejections(diff_phos, alpha = 10^-0.05, lfc = 1)
    return(sig_phos)
}

#' @export
dep2_rnaseq <- function(df, tbl_name, rv) {
    unique_dds <- df
    ecols <- match(trimws(rv$data_cols[[tbl_name]]), colnames(unique_dds))
    data_mat <- as.matrix(unique_dds[, ecols])
    rownames(data_mat) <- paste0(unique_dds$Gene_Name, "_", unique_dds$Gene_ID)
    dds <- DEP2::make_dds_parse(data_mat, mode = "delim", sep = "_")
    dds <- DEP2::filter_se(dds, fraction = 0.3, thr = 1)
    diff <- DEP2::test_diff_deg(dds, type = "all")
    sig <- DEP2::add_rejections(diff, alpha = 10^-0.05, lfc = 1)
    return(sig)
}
