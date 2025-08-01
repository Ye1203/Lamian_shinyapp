color_selected <- function(color_length) {
  color_total <- c(
    "#e6194b","#ffe119","#46f0f0","#f58231","#bcf60c", 
    "#ff00ff","#9a6324","#fffac8","#e6beff","#00bfff", 
    "#ffd8b1","#00ff7f","#f5a9bc","#1e90ff","#ffa500",
    "#98fb98","#911eb4","#afeeee","#fa8072","#9acd32",
    "#3cb44b","#000075","#808000","#cd5c5c","#dda0dd",
    "#40e0d0","#ff69b4","#8a2be2","#c71585","#5f9ea0",
    "#dc143c","#87cefa","#ff6347","#9932cc","#00ced1",
    "#ff4500","#6a5acd","#b0e0e6","#d2691e","#a9a9f5",
    "#adff2f","#8b0000","#7fffd4","#00fa9a","#ba55d3",
    "#2e8b57","#ffdab9","#b22222","#ffe4e1","#7b68ee"
  )
  
  if (color_length <= length(color_total)) {
    return(color_total[1:color_length])
  }
  
  warning("groups is larger than 50, color will randomly select")
  palette_fn <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
  extra_colors <- palette_fn(color_length - length(color_total))
  colors <- c(color_total, extra_colors)
  return(colors)
}

draw_heatmap <- function(xde_result, gene_list, cluster_method, scale_method, low_color, mid_color, high_color, cluster_number) {
  warning_msg <- NA
  set.seed(123)
  valid_genes <- intersect(gene_list, rownames(xde_result$expr))
  if(length(valid_genes) == 0) {
    stop("None of the input genes exist in the dataset")
  }
  
  missing_genes <- setdiff(gene_list, valid_genes)
  if(length(missing_genes) > 0) {
    warning_msg <- paste(
      length(missing_genes), 
      "genes not found:",
      paste(head(missing_genes, 5), collapse = ", "),
      ifelse(length(missing_genes) > 5, "...", "")
    )
  }
  
  col_fun <- circlize::colorRamp2(
    c(-2, 0, 2), 
    c(low_color, mid_color, high_color) 
  )
  
  # Get raw data
  raw0 <- xde_result$populationFit[[1]][valid_genes, , drop = FALSE]
  raw1 <- xde_result$populationFit[[2]][valid_genes, , drop = FALSE]
  
  if (scale_method == "Together") {
    raw_both <- cbind(raw0, raw1)
    fit_both_scaled <- t(scale(t(raw_both))) 
    ncol_raw0 <- ncol(raw0)
    scale0 <- fit_both_scaled[, 1:ncol_raw0, drop = FALSE]
    scale1 <- fit_both_scaled[, (ncol_raw0 + 1):ncol(fit_both_scaled), drop = FALSE]
  } else {
    scale0 <- t(scale(t(raw0)))
    scale1 <- t(scale(t(raw1)))
  }
  
  # Cluster the data
  if(cluster_method == "kmeans"){
    km_res <- kmeans(fit_both_scaled, centers = cluster_number, iter.max = 1000)
    ord <- order(km_res$cluster)
    clusters <- km_res$cluster[ord]
  } else {
    dist_mat <- dist(fit_both_scaled)
    hc_res <- hclust(dist_mat, method = "complete")
    clusters <- cutree(hc_res, k = cluster_number)
    ord <- hc_res$order
  }
  
  # Reorder all data
  gene_names <- valid_genes[ord]
  scale0 <- scale0[ord, , drop = FALSE]
  scale1 <- scale1[ord, , drop = FALSE]
  raw0 <- raw0[ord, , drop = FALSE]
  raw1 <- raw1[ord, , drop = FALSE]
  clusters <- clusters[ord]
  
  # Create the output data frame with requested column names
  df <- data.frame(
    GENE = gene_names,
    Cluster_Number = clusters,
    stringsAsFactors = FALSE
  )
  
  # Add scaled data with scale0_ prefix
  scale0_df <- as.data.frame(scale0)
  colnames(scale0_df) <- paste0("scale0_", seq_len(ncol(scale0)))
  df <- cbind(df, scale0_df)
  
  # Add scaled data with scale1_ prefix
  scale1_df <- as.data.frame(scale1)
  colnames(scale1_df) <- paste0("scale1_", seq_len(ncol(scale1)))
  df <- cbind(df, scale1_df)
  
  # Add raw data with raw0_ prefix
  raw0_df <- as.data.frame(raw0)
  colnames(raw0_df) <- paste0("raw0_", seq_len(ncol(raw0)))
  df <- cbind(df, raw0_df)
  
  # Add raw data with raw1_ prefix
  raw1_df <- as.data.frame(raw1)
  colnames(raw1_df) <- paste0("raw1_", seq_len(ncol(raw1)))
  df <- cbind(df, raw1_df)
  
  # Create cluster colors
  cluster_colors <- setNames(
    color_selected(length(unique(clusters))),
    sort(unique(clusters))
  )
  
  # Create heatmap
  row_ha <- rowAnnotation(
    Cluster = as.factor(clusters),
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE
  )
  
  ht0 <- Heatmap(
    scale0,
    name = "male",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    col = col_fun,
    left_annotation = row_ha
  )
  
  ht1 <- Heatmap(
    scale1,
    name = "female",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    col = col_fun
  )
  
  p <- draw(ht0 + ht1, merge_legend = TRUE, heatmap_legend_side = "bottom")
  
  list(
    plot = p,
    df = df,
    warning_message = warning_msg
  )
}