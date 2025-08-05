
step4 <- function(obj, sce, non_zero_num, lower_quantile, upper_quantile, selected_lineage){
  expr <- GetAssayData(obj, layer = "data")
  
  originial_cells_num <- ncol(expr)
  expr_count <- rowSums(expr > 0)
  gene_list <- names(expr_count)[expr_count >= non_zero_num]
  
  pseudotime <- slingPseudotime(sce, na = TRUE)[,selected_lineage]
  pseudotime <- pseudotime[!is.na(pseudotime)]
  names(pseudotime) <- colnames(sce)
  cell_list <- intersect(colnames(obj),names(pseudotime)[!is.na(pseudotime)])
  
  pseudotime <- pseudotime[cell_list]
  Q1 <- quantile(pseudotime, lower_quantile/100, na.rm = TRUE)  # swapped upper/lower quantile
  Q3 <- quantile(pseudotime, upper_quantile/100, na.rm = TRUE)
  inlier <- pseudotime >= Q1 & pseudotime <= Q3
  
  keep_cells <- names(pseudotime)[inlier]
  sce_sub <- sce[, keep_cells]
  obj_sub <- subset(obj, cells = keep_cells, features = gene_list)
  return(list(
    obj = obj_sub,
    sce_sub = sce_sub,
    non_zero_num = non_zero_num,
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile
  ))
}