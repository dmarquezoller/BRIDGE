#' @export
dep2_proteomics <- function(df, tbl_name, rv) {
  unique_pg <- DEP2::make_unique(df, names = "Gene_Name", ids = "Protein_ID", delim = ";")
  ecols <- match(trimws(rv$time_cols[[tbl_name]]), colnames(unique_pg))
  se_pg <- DEP2::make_se_parse(unique_pg, columns = ecols, mode = "delim", sep = "_", remove_prefix = T, log2transform = T)
  filter_pg <- DEP2::filter_se(se_pg, thr = 1, fraction = 0.3)
  impute_pg <- DEP2::impute(filter_pg, fun = "MinDet")
  norm_pg <- DEP2::normalize_vsn(impute_pg)
  diff_pg <- DEP2::test_diff(norm_pg, type = "all", fdr.type = "BH")
  sig_pg <- DEP2::add_rejections(diff_pg, alpha = 10^-0.05, lfc = 1)
  return(sig_pg)
}

#' @export
dep2_phosphoproteomics <- function(df, tbl_name, rv) {
  df <- df[!grepl("p0", df$pepG), ] # Remove the p0
  unique_phos <-DEP2::make_unique(df, names = "pepG", ids = "Protein_ID", delim = ";")
  ecols <- match(trimws(rv$time_cols[[tbl_name]]), colnames(unique_phos))
  se_phos <- DEP2::make_se_parse(unique_phos, columns = ecols, mode = "delim", sep = "_", remove_prefix = T, log2transform = T)
  filter_phos <- DEP2::filter_se(se_phos)
  impute_phos <- DEP2::impute(filter_phos, fun = "MinDet")
  norm_phos <- DEP2::normalize_vsn(impute_phos)
  diff_phos <- DEP2::test_diff(norm_phos, type = "all", fdr.type = "BH")
  sig_phos <- DEP2::add_rejections(diff_phos, alpha = 10^-0.05, lfc = 1)
  return(diff_phos)
}

#' @export
dep2_rnaseq <- function(df, tbl_name, rv) {
  unique_dds <- df
  ecols <- match(trimws(rv$time_cols[[tbl_name]]), colnames(unique_dds))
  data_mat <- as.matrix(unique_dds[, ecols])
  rownames(data_mat) <- unique_dds$Gene_ID
  dds <- DEP2::make_dds_parse(data_mat, mode = "delim", sep = "_")
  dds <- DEP2::filter_se(dds, fraction = 0.3, thr = 1)
  diff <- DEP2::test_diff_deg(dds, type = 'all')
  sig <- DEP2::add_rejections(diff, alpha = 10^-0.05, lfc = 1)
  return(sig)
}