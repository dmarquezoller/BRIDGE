#' @export
plot_custom_tc_cluster <- function(ht_mat, exp_design, cluster_vec,
                                   groupby = "condition", group_order = NULL,
                                   col_limit = 4, row_font_size = 5, col_font_size = 5,
                                   heatmap_width = 3, heatmap_height = 5,
                                   color = "RdBu") {
  
  # Preprocess inputs
  color <- match.arg(color)
  if (is.null(group_order)) {
    group_order <- unique(exp_design[[groupby]])
  }
  
  # Ensure correct matrix order
  ht_mat <- ht_mat[, which(exp_design[[groupby]] %in% group_order)]
  ht_mat <- ht_mat - rowMeans(ht_mat)
  exp_design2 <- exp_design[which(exp_design[[groupby]] %in% group_order), ]
  
  # Compute tcadata matrix for line plots
  ht_mat2 <- group_order %>% map_df(function(g) {
    rowMeans(ht_mat[, exp_design2[[groupby]] == g, drop = FALSE])
  }) %>% t()
  colnames(ht_mat2) <- group_order
  rownames(ht_mat2) <- rownames(ht_mat)
  
  # Ensure cluster_vec is a factor, ordered by ht_mat rows
  cluster_vec <- factor(unlist(cluster_vec))
  cluster_vec <- cluster_vec[rownames(ht_mat)]
  
  # Dummy membership values for coloring
  k <- length(unique(cluster_vec))
  set.seed(1234)
  mebership <- runif(nrow(ht_mat), 0.1, 0.9)
  names(mebership) <- rownames(ht_mat)
  
  # Color scales
  col_fun_line <- circlize::colorRamp2(c(0.01, 0.5, 0.9), c("#54c556", "#d3ca48", "#cd3634"))
  col_fun_heatmap <- circlize::colorRamp2(seq(-col_limit, col_limit, length.out = 11),
                                          rev(RColorBrewer::brewer.pal(11, color)))
  
  ngroup <- length(group_order)
  
  # Line plot function
  panel_fun2 <- function(index, nm) {
    pushViewport(viewport(xscale = c(0, (2 * ngroup + 1)), yscale = c(-2, 2)))
   grid::grid.rect()
    
    gby <- ComplexHeatmap::annotation_axis_grob(
      at = c(-1.7, 0, 1.7),
      labels = c(-1.7, 0, 1.7),
      labels_rot = 0,
      side = "left", facing = "inside",
      gp =ggfun::gpar(fontsize = 5)
    )
    
   grid::grid.polyline(
      y = as.vector(t((ht_mat2[index, , drop = FALSE] / 4) + 0.5)),
      x = rep(seq(from = 2 / (2 * ngroup + 2), to = (2 * ngroup) / (2 * ngroup + 1), length.out = ngroup), length(index)),
      id = rep(1:length(index), each = ngroup),
      gp =ggfun::gpar(col = col_fun_line(mebership[index]), lwd = 1, alpha = 0.5)
    )
    
    if (cluster_vec[index[1]] == levels(cluster_vec)[k]) {
      gbx <- ComplexHeatmap::annotation_axis_grob(
        at = seq(2, 2 * ngroup, by = 2),
        labels = colnames(ht_mat2),
        labels_rot = 45,
        side = "bottom", facing = "outside",
        gp =ggfun::gpar(fontsize = col_font_size)
      )
     grid::grid.draw(gbx)
    }
    
   grid::grid.draw(gby)
    popViewport()
  }
  
  # Row annotation with trend lines
  anno <- ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_link(
    align_to = cluster_vec,
    which = "row",
    panel_fun = panel_fun2,
    side = "right",
    size = unit((heatmap_height * 5 / k) - 0.25, "cm"),
    gap = unit(0, "cm"),
    width = unit(heatmap_width * 1.7, "cm")
  ))
  
  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    ht_mat,
    heatmap_width = unit(heatmap_width * 5, "cm"),
    heatmap_height = unit(heatmap_height * 5, "cm"),
    col = col_fun_heatmap,
    split = cluster_vec,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_dend_width = unit(10, "mm"),
    cluster_row_slices = TRUE,
    cluster_column_slices = FALSE,
    row_gap = unit(1, "mm"),
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = "Centered value", title_position = "leftcenter-rot"),
    row_names_gp =ggfun::gpar(fontsize = row_font_size),
    column_names_gp =ggfun::gpar(fontsize = col_font_size),
    right_annotation = anno
  )
  
  return(ht)
}