
step3 <- function(obj, sample, clusters, reduction, start_cluster) {
  expr_mat <- GetAssayData(obj, layer = "data")
  sce <- SingleCellExperiment(
    assays = list(logcounts = expr_mat))
  sce$sample <- obj@meta.data[[sample]]
  sce$cluster <- obj@meta.data[[clusters]]
  reducedDims(sce)$Reduction <- Embeddings(obj, reduction)
  reducedDims(sce)$UMAP <- Embeddings(obj, "umap")
  sce <- slingshot(
    sce,
    clusterLabels = "cluster",
    reducedDim    = "Reduction",
    start.clus    = start_cluster  
  )
  
  # # plot1
  # xvar <- paste0(reduction, "_1")
  # yvar <- paste0(reduction, "_2")
  # 
  # rd_df <- as.data.frame(reducedDims(sce)$Reduction[, 1:2])
  # colnames(rd_df) <- c(xvar, yvar)
  # rd_df$cluster <- sce$cluster
  # 
  # lineages <- slingCurves(sce)
  # curve_data <- data.frame()
  # for (i in seq_along(lineages)) {
  #   curve_i <- as.data.frame(lineages[[i]]$s[,1:2])
  #   colnames(curve_i) <- c(xvar, yvar)
  #   curve_i$lineage <- i
  #   curve_data <- rbind(curve_data, curve_i)
  # }
  # curve_data_red <- subset(curve_data, lineage == 1)
  # curve_data_black <- subset(curve_data, lineage != 1)
  # 
  # umap_plot <- ggplot(rd_df, aes(x = !!sym(xvar), y = !!sym(yvar), color = cluster)) +
  #   geom_point(size = 0.5) +
  #   labs(x = xvar, y = yvar, color = "Cluster Labels") +  
  #   theme_minimal() +
  #   theme(
  #     legend.title = element_text(size = 0),  
  #     legend.text = element_text(size = 8),  
  #     plot.title = element_blank(),  
  #     axis.title = element_text(size = 10, face = "bold"),  
  #     axis.text = element_text(size = 8)  
  #   )
  # 
  # plot1 <- umap_plot + 
  #   geom_path(data = curve_data_black,
  #             aes(x = !!sym(xvar), y = !!sym(yvar), group = lineage),
  #             color = "black",
  #             size = 1) +
  #   geom_path(data = curve_data_red,
  #             aes(x = !!sym(xvar), y = !!sym(yvar), group = lineage),
  #             color = "red",
  #             size = 1.5)
  # 
  # # plot 2
  # df <- data.frame(
  #   Dim_1     = reducedDims(sce)$Reduction[,1],
  #   Dim_2     = reducedDims(sce)$Reduction[,2],
  #   cluster    = colData(sce)$cluster,
  #   pseudotime = colData(sce)$slingPseudotime_1
  # )
  # 
  # plot2 <- ggplot(df, aes(Dim_1, Dim_2, color = pseudotime)) +
  #   geom_point(size = 1.5) +
  #   scale_color_viridis_c(option = "D") +
  #   labs(title = "Slingshot colored by pseudotime", color = "pseudotime") +
  #   theme_classic()
  
  return(sce)
}