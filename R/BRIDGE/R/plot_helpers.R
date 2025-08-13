# open a temp PNG device, evaluate expr (which should DRAW), close, return file path + dims
.with_png_device <- function(expr, width = 1200, height = 900, res = 144) {
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf, width = width, height = height, res = res)
  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    try(graphics::par(op), silent = TRUE)
    try(grDevices::dev.off(), silent = TRUE)
  }, add = TRUE)
  force(expr)
  list(src = tf, width = width, height = height)
}

# 1) Run NbClust without leaking devices; return the chosen k
sandbox_nbclust_k <- function(mat, min.nc = 2, max.nc = 10) {
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf, width = 1200, height = 900, res = 144)
  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    try(graphics::par(op), silent = TRUE)
    try(grDevices::dev.off(), silent = TRUE)
    try(unlink(tf), silent = TRUE)
  }, add = TRUE)
  elbow <- NbClust::NbClust(as.matrix(mat),
                            distance = "euclidean",
                            min.nc   = min.nc,
                            max.nc   = max.nc,
                            method   = "kmeans")
  as.numeric(names(sort(table(elbow$Best.nc[1, ]), decreasing = TRUE)[1]))
}

# 2) Draw a DEP heatmap to a PNG and return list(src, width, height)
draw_heatmap_png <- function(dep_out_cut, clustering_enabled, k, width = 1200, height = 900, res = 144) {
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf, width = width, height = height, res = res)
  on.exit({ try(grDevices::dev.off(), silent = TRUE) }, add = TRUE)
  if (isTRUE(clustering_enabled) && is.finite(k) && !is.na(k)) {
    DEP2::plot_heatmap(dep_out_cut, kmeans = TRUE, k = k)
  } else {
    DEP2::plot_heatmap(dep_out_cut)
  }
  list(src = tf, width = width, height = height)
}
