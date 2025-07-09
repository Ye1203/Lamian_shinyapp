
step3 <- function(obj, sample, clusters, reduction, start_cluster) {
  expr_mat <- GetAssayData(obj, slot = "data")
  sce <<- SingleCellExperiment(
    assays = list(logcounts = expr_mat))
  sce$sample <- obj@meta.data[[sample]]
  sce$cluster <- obj@meta.data[[clusters]]
  reducedDims(sce)$Reduction <- Embeddings(obj, reduction)
  
  sce <- slingshot(
    sce,
    clusterLabels = "cluster",
    reducedDim    = "Reduction",
    start.clus    = start_cluster  
  )
  
  df <- data.frame(
    Dim_1     = reducedDims(sce)$Reduction[,1],
    Dim_2     = reducedDims(sce)$Reduction[,2],
    cluster    = colData(sce)$cluster,
    pseudotime = colData(sce)$slingPseudotime_1
  )
  plot1 <- ggplot(df, aes(Dim_1, Dim_2, color = cluster)) +
    geom_point(size = 1.5) +
    labs(title = "Slingshot colored by cluster") +
    theme_classic()
  
  plot2 <- ggplot(df, aes(Dim_1, Dim_2, color = pseudotime)) +
    geom_point(size = 1.5) +
    scale_color_viridis_c(option = "D") +
    labs(title = "Slingshot colored by pseudotime", color = "pseudotime") +
    theme_classic()
  
  return(list(
    sce = sce,
    plot1 = plot1,
    plot2 = plot2))
}