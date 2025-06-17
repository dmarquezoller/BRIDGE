#' @export
total_intensity_normalization <- function(raw_data, pseudocount = 1) {
  # Identify sample columns (numeric only)
  sample_cols <- sapply(raw_data, is.numeric)
  
  # Extract just the expression matrix
  expr_data <- raw_data[, sample_cols]
  
  # Step 1: Total intensity per sample
  total_intensity <- colSums(expr_data, na.rm = TRUE)
  
  # Step 2: Max total intensity
  max_total <- max(total_intensity)
  
  # Step 3: Normalize
  norm_data <- sweep(expr_data, 2, total_intensity, FUN = "/") * max_total
  
  # Step 4: Log2 transform with pseudocount
  norm_data_log2 <- log2(norm_data + pseudocount)
  
  # Step 5: Combine with metadata columns
  meta_cols <- raw_data[, !sample_cols]
  processed_data <- cbind(meta_cols, norm_data_log2)
  
  return(processed_data)
}

#' @export
median_normalization <- function(raw_data, pseudocount = 1) {
  # Identify sample columns (numeric only)
  sample_cols <- sapply(raw_data, is.numeric)
  
  # Extract numeric expression matrix
  expr_data <- raw_data[, sample_cols]
  
  # Step 1: Median intensity per sample
  sample_medians <- apply(expr_data, 2, median, na.rm = TRUE)
  
  # Step 2: Maximum median across samples
  max_median <- max(sample_medians)
  
  # Step 3: Normalize by sample median and scale to max median
  norm_data <- sweep(expr_data, 2, sample_medians, FUN = "/") * max_median
  
  # Step 4: Log2 transform with pseudocount
  norm_data_log2 <- log2(norm_data + pseudocount)
  
  # Step 5: Combine with metadata columns
  meta_cols <- raw_data[, !sample_cols]
  processed_data <- cbind(meta_cols, norm_data_log2)
  
  return(processed_data)
}

#' @export
log2_transform <- function(raw_data, pseudocount = 1) {
  # Identify sample columns (numeric only)
  sample_cols <- sapply(raw_data, is.numeric)
  
  # Extract and transform expression data
  expr_data <- raw_data[, sample_cols]
  log2_data <- log2(expr_data + pseudocount)
  
  # Combine with metadata
  meta_cols <- raw_data[, !sample_cols]
  processed_data <- cbind(meta_cols, log2_data)
  
  return(processed_data)
}

#' @export
fpkm_normalization <- function(raw_data, annotation) {
  # Compute gene length
  annotation$Gene_Length <- annotation$Gene_End - annotation$Gene_Start + 1

  # Standardize Gene_Name to lowercase for matching
  raw_data$Gene_Name_lower <- tolower(raw_data$Gene_Name)
  annotation$Gene_Name_lower <- tolower(annotation$Gene_Name)

  # Match by lowercase gene name
  matched_idx <- match(raw_data$Gene_Name_lower, annotation$Gene_Name_lower)
  gene_lengths <- annotation$Gene_Length[matched_idx]

  # Check for missing values
  if (any(is.na(gene_lengths))) {
    warning("Some genes could not be matched and will be excluded.")
  }

  # Filter out unmatched rows
  keep <- !is.na(gene_lengths)
  raw_data <- raw_data[keep, ]
  gene_lengths <- gene_lengths[keep]

  # Extract expression matrix
  sample_cols <- sapply(raw_data, is.numeric)
  expr_data <- as.matrix(raw_data[, sample_cols])

  # Calculate RPKM
  rpkm_data <- edgeR::rpkm(expr_data, gene.length = gene_lengths)

  # Combine with metadata
  meta_cols <- raw_data[, !sample_cols & !names(raw_data) %in% "Gene_Name_lower"]
  processed_data <- cbind(meta_cols, as.data.frame(rpkm_data))

  return(processed_data)
}

#' @export
tmm_normalization <- function(raw_data) {
  # Identify numeric columns (assumed to be raw counts)
  sample_cols <- sapply(raw_data, is.numeric)
  expr_data <- as.matrix(raw_data[, sample_cols])

  # Create DGEList object
  dge <- edgeR::DGEList(counts = expr_data)

  # Compute TMM normalization factors
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  # Get library sizes and normalization factors
  lib_sizes <- dge$samples$lib.size
  norm_factors <- dge$samples$norm.factors
  effective_lib_sizes <- lib_sizes * norm_factors

  # Scale raw counts by effective library size (TMM-normalized counts)
  scaling_factors <- effective_lib_sizes / mean(effective_lib_sizes)
  norm_counts <- t(t(expr_data) / scaling_factors)

  # Combine with metadata
  meta_cols <- raw_data[, !sample_cols]
  processed_data <- cbind(meta_cols, as.data.frame(norm_counts))

  return(processed_data)
}

#' @export
cpm_normalization <- function(raw_data) {
  # Identify numeric columns (assumed to be expression data)
  sample_cols <- sapply(raw_data, is.numeric)
  expr_data <- as.matrix(raw_data[, sample_cols])

  # Create DGEList object
  dge <- edgeR::DGEList(counts = expr_data)

  # Calculate TMM normalization factors
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  # Get normalized CPM values (log2 optional, set log = FALSE for raw CPMs)
  norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

  # Combine with metadata (non-numeric columns)
  meta_cols <- raw_data[, !sample_cols]
  processed_data <- cbind(meta_cols, as.data.frame(norm_counts))

  return(processed_data)
}

#' @export
tpm_normalization <- function(fpkm_data) {
  # Identify numeric columns (assumed to be FPKM values)
  sample_cols <- sapply(fpkm_data, is.numeric)
  fpkm_matrix <- as.matrix(fpkm_data[, sample_cols])

  # Apply the TPM conversion column-wise
  tpm_matrix <- apply(fpkm_matrix, 2, function(x) {
    exp(log(x) - log(sum(x, na.rm = TRUE)) + log(1e6))
  })

  # Combine with metadata
  meta_cols <- fpkm_data[, !sample_cols]
  processed_data <- cbind(meta_cols, as.data.frame(tpm_matrix))

  return(processed_data)
}


