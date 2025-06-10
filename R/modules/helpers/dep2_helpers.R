library(DEP2)

dep2_proteomics <- function(df, tbl_name, rv) {
  unique_pg <- make_unique(df, names = "Gene_Name", ids = "Protein_ID", delim = ";")
  ecols <- match(trimws(rv$time_cols[[tbl_name]]), colnames(unique_pg))
  se_pg <- DEP2::make_se_parse(unique_pg, columns = ecols, mode = "char", chars = 1, remove_prefix = T, log2transform = T)
  filter_pg <- filter_se(se_pg, thr = 1, fraction = 0.3)
  norm_pg <- normalize_vsn(filter_pg)
  diff_pg <- test_diff(norm_pg, type = "all", fdr.type = "BH")
  return(diff_pg)
}

dep2_phosphoproteomics <- function(df, tbl_name, rv) {
  df <- df[!grepl("p0", df$pepG), ] # Remove the p0
  unique_phos <- make_unique(df, names = "pepG", ids = "Protein_ID", delim = ";")
  ecols <- match(trimws(rv$time_cols[[tbl_name]]), colnames(unique_phos))
  se_phos <- DEP2::make_se_parse(unique_phos, columns = ecols, mode = "char", chars = 1, remove_prefix = T, log2transform = T)
  filter_phos <- filter_se(se_phos)
  norm_phos <- normalize_vsn(filter_phos)
  diff_phos <- test_diff(norm_phos, type = "all", fdr.type = "BH")
  return(diff_phos)
}

dep2_rnaseq <- function(df, tbl_name, rv) {
  unique_dds <- df
  ecols <- match(trimws(rv$time_cols[[tbl_name]]), colnames(unique_dds))
  data_mat <- as.matrix(unique_dds[, ecols])
  rownames(data_mat) <- unique_dds$Gene_ID
  dds <- make_dds_parse(data_mat, mode = 'char', chars = 1)
  dds <- filter_se(dds, fraction = 0.3, thr = 1)
  diff <- test_diff_deg(dds, type = 'all')
  return(diff)
}