#testing step1(data = obj, is_harmony, group_by_vars, dims_use)
step1 <- function(data,
                  is_harmony,
                  group_by_vars,
                  dims_use){
  DefaultAssay(data) <- "RNA"
  data <- NormalizeData(data, normalization.method = "LogNormalize", verbose = FALSE)
  if(is_harmony){
    data <- RunHarmony(
      object = data,
      group.by.vars = group_by_vars,
      reduction.use = "pca",
      dims.use = 1:dims_use,
      assay.use = "RNA",
      reduction.save = "harmony",
      project.dim = TRUE,
      verbose = FALSE
    )
  }
  return(data)
}