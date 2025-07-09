
step4 <- function(obj, sce, non_zero_num, lower_quantile, upper_quantile, IQR_coefficient){
  obj <<- obj
  sce <<- sce
  expr <- GetAssayData(obj, layer = "data")
  originial_cells_num <- ncol(expr)
  expr_count <- rowSums(expr > 0)
  gene_list <- names(expr_count)[expr_count >= non_zero_num]
  
  pseudotime <- slingPseudotime(sce, na = TRUE)[,1]
  names(pseudotime) <- colnames(sce)
  cell_list <- intersect(colnames(obj),names(pseudotime)[!is.na(pseudotime)])
  
  pseudotime <- pseudotime[cell_list]
  Q1 <- quantile(pseudotime, lower_quantile, na.rm = TRUE)  # swapped upper/lower quantile
  Q3 <- quantile(pseudotime, upper_quantile, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  inlier <- pseudotime >= (Q1 - IQR_coefficient * IQR_val) & pseudotime <= (Q3 + IQR_coefficient * IQR_val)
  
  keep_cells <- names(pseudotime)[inlier]
  obj_sub <- subset(obj, cells = keep_cells, features = gene_list)
  return(list(
    obj = obj_sub,
    Q1 = Q1,
    Q3 = Q3,
    IQR = IQR_val
  ))
}