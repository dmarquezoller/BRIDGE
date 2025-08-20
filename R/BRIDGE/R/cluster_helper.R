safe_nbclust <- function(mat, k_min = 2, k_max = 10,
                         index_primary = "ch", index_fallback = "silhouette",
                         seed = 123L) {
  mat <- as.matrix(mat)

  # Replace NaN/Inf with NA, then drop bad rows/cols and zero-variance ones
  mat[!is.finite(mat)] <- NA_real_
  row_ok <- rowSums(is.finite(mat)) >= 2 &
            apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  col_ok <- colSums(is.finite(mat)) >= 2 &
            apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
  mat <- mat[row_ok, col_ok, drop = FALSE]

  # Guard: need at least 2x2 and k_max ≤ nrow - 1
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_integer_)
  k_max <- max(k_min, min(k_max, nrow(mat) - 1))

  # Avoid device warnings inside futures
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)

  # Try primary index ("ch"), then fallback ("silhouette")
  try_indices <- c(index_primary, index_fallback)
  for (idx in unique(try_indices)) {
    set.seed(seed)
    res <- try(
      NbClust::NbClust(mat, distance = "euclidean",
                       min.nc = k_min, max.nc = k_max,
                       method  = "kmeans", index = idx),
      silent = TRUE
    )
    if (!inherits(res, "try-error")) {
      best <- res$Best.nc
      if (is.matrix(best)) {
        kv <- as.integer(best[1, ])
        return(as.integer(names(which.max(table(kv)))))
      } else {
        return(as.integer(best[1]))
      }
    }
  }

  # Last-resort: simple elbow on total within-ss (robust and fast)
  set.seed(seed)
  ks  <- k_min:k_max
  wss <- vapply(ks, function(k) {
    stats::kmeans(mat, centers = k, nstart = 5, iter.max = 50)$tot.withinss
  }, numeric(1))
  # Choose elbow by maximum second difference
  if (length(wss) >= 3) {
    dd <- diff(diff(wss))
    return(ks[which.max(c(NA, dd, NA))])
  }
  ks[which.min(wss)]
}

safe_row_scale <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sd <- apply(m, 1, function(x) stats::sd(x, na.rm = TRUE))
  sd[!is.finite(sd) | sd == 0] <- 1
  m <- sweep(m, 1, mu, FUN = "-")
  sweep(m, 1, sd, FUN = "/")
}

