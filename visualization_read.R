visualization_read <- function(result_folder){
  xde_result <- readRDS(file.path(result_folder,"xde_result.RDS"))
  sce <- readRDS(file.path(result_folder,"sce.RDS"))
  expr <- readRDS(file.path(result_folder,"expr.RDS"))
  xde_result$expr <- expr
  stat <- xde_result$statistics
  stat <- as.data.frame(stat)
  stat <- stat[order(stat$fdr.overall), ]
  xde_result$populationFit <-
    getPopulationFit(xde_result, gene = rownames(stat), type = 'variable')
  xde_result$covariateGroupDiff <-
    getCovariateGroupDiff(testobj = xde_result, gene = rownames(stat))
  return(list(xde_result = xde_result,
              stat = stat,
              sce = sce))

  }