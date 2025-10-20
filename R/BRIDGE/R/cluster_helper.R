safe_nbclust <- function(mat, k_min = 2, k_max = 10,
                         index_primary = "ch", index_fallback = "silhouette",
                         seed = 123L) {
    mat <- as.matrix(mat)
    # sanitize: finite & non-constant rows/cols
    mat[!is.finite(mat)] <- NA_real_
    row_ok <- rowSums(is.finite(mat)) >= 2 &
        apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
    col_ok <- colSums(is.finite(mat)) >= 2 &
        apply(mat, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
    mat <- mat[row_ok, col_ok, drop = FALSE]

    if (nrow(mat) < 2 || ncol(mat) < 2) {
        return(NA_integer_)
    }

    # max allowed k is nrow-1
    k_max_eff <- min(k_max, nrow(mat) - 1)
    if (k_max_eff < k_min) {
        return(NA_integer_)
    }

    # avoid worker device warnings
    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp)
    on.exit(
        {
            dev.off()
            unlink(tmp)
        },
        add = TRUE
    )

    # try the NbClust indices first
    for (idx in unique(c(index_primary, index_fallback))) {
        set.seed(seed)
        res <- try(
            NbClust::NbClust(mat,
                distance = "euclidean",
                min.nc = k_min, max.nc = k_max_eff,
                method = "kmeans", index = idx
            ),
            silent = TRUE
        )
        if (!inherits(res, "try-error")) {
            best <- res$Best.nc
            return(if (is.matrix(best)) {
                kv <- as.integer(best[1, ])
                as.integer(names(which.max(table(kv))))
            } else {
                as.integer(best[1])
            })
        }
    }

    # last resort: elbow from kmeans() with per-k clamps & tryCatch
    ks <- seq.int(k_min, k_max_eff)
    ks <- ks[ks >= 1 & ks < nrow(mat)]
    if (!length(ks)) {
        return(NA_integer_)
    }

    set.seed(seed)
    wss <- vapply(ks, function(k) {
        out <- try(stats::kmeans(mat, centers = k, nstart = 5, iter.max = 50),
            silent = TRUE
        )
        if (inherits(out, "try-error")) NA_real_ else out$tot.withinss
    }, numeric(1))

    if (!any(is.finite(wss))) {
        return(NA_integer_)
    }

    # choose elbow by max second difference, else min WSS
    if (length(wss) >= 3) {
        dd <- diff(diff(wss))
        # align index; pad NAs at ends
        pick <- which.max(dd)
        if (length(pick) && is.finite(dd[pick])) {
            return(ks[pick + 1L])
        }
    }
    ks[which.min(wss)]
}

# robust row scaling: never produces NaN rows
safe_row_scale <- function(m) {
    m <- as.matrix(m)
    mu <- matrixStats::rowMeans2(m, na.rm = TRUE)
    sdv <- matrixStats::rowSds(m, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] <- 1 # prevent divide-by-zero
    m <- sweep(m, 1, mu, FUN = "-")
    sweep(m, 1, sdv, FUN = "/")
}
