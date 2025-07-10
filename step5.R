step5 <- function(data, sce, output_file_path, design, maximumknotallowed, permuiter, project_name, email_address){
  if (!startsWith(output_file_path, "/")) {
    output_file_path <- paste0("/", output_file_path)
  }
  
  if (!dir.exists(output_file_path)) {
    dir.create(output_file_path, recursive = TRUE)
  }
  
  data@assays$SCT <- NULL
  pseudotime <- slingPseudotime(sce, na = TRUE)[,1]
  names(pseudotime) <- colnames(sce)
  expr <- as.matrix(GetAssayData(data, layer = "data"))
  cellanno <- data.frame(
    Cell = colnames(expr),
    Sample  = as.character(data@meta.data$sample_id),
    stringsAsFactors = FALSE
  )
  pseudotime <- pseudotime[names(pseudotime) %in% cellanno$Cell]
  
  rownames(design) <- design$Sample
  design$Sample <- NULL
  
  saveRDS(expr, file.path(output_file_path,"expr.RDS"))
  saveRDS(pseudotime, file.path(output_file_path, "pseudotime.RDS"))
  saveRDS(design, file.path(output_file_path, "design.RDS"))
  saveRDS(cellanno, file.path(output_file_path, "cellanno.RDS"))
  
  r_script_path <- file.path(output_file_path, "lamian.R")
  
  lamian_r_code <- sprintf(
    'options(warn=-1) # suppress warning messages
library(devtools)
devtools::load_all("/projectnb/wax-es/00_shinyapp/Lamian/renv/library/linux-almalinux-8.10/R-4.4/x86_64-pc-linux-gnu/Lamian")
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(slingshot)
library(SingleCellExperiment)
source("/projectnb/wax-es/00_shinyapp/Lamian/renv/activate.R")

send_email <- function(to, subject, body) {
  f <- tempfile()
  writeLines(body, f)
  cmd <- sprintf("mail -s \\"%%s\\" %%s < %%s", subject, to, f)
  system(cmd)
  unlink(f)
}

output_file_path <- "%s"
email_address <- "%s"
maximumknotallowed <- %d
permuiter <- %d
result_folder <- output_file_path

cellanno <- readRDS(file.path(output_file_path,"cellanno.RDS"))
expr <- readRDS(file.path(output_file_path,"expr.RDS"))
design <- readRDS(file.path(output_file_path,"design.RDS"))
pseudotime <- readRDS(file.path(output_file_path,"pseudotime.RDS"))

start_time <- Sys.time()

res <- tryCatch({
  xde_result <- lamian_test(
    expr = expr,
    cellanno = cellanno,
    pseudotime = pseudotime,
    design = design,
    test.type = "variable",
    maxknotallowed = maximumknotallowed,
    testvar = 2,
    permuiter = permuiter,
    test.method = "permutation",
    ncores = 28,
    verbose.output = TRUE
  )
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  saveRDS(xde_result, file.path(result_folder,"xde_result.RDS"))
  
  body <- c(
    "Hi,",
    "",
    sprintf("The procedure has been completed successfully, the result is under folder %%s. Spending %%s", output_file_path, elapsed_time),
    "",
    "Best,",
    "Bingtian"
  )
  send_email(subject = "Program Completed", body = body, to = email_address)
  
  list(success = TRUE, elapsed = elapsed_time)
}, error = function(e) {
  err_msg <- conditionMessage(e)
  body <- c(
    "Hi,",
    "",
    sprintf("There is an error (%%s). You can copy this email and send to btye@bu.edu.", err_msg),
    "",
    sprintf("output_file_path = %%s", output_file_path),
    "",
    "Best,",
    "Bingtian"
  )
  send_email(subject = "Program Completed", body = body, to = email_address)
  list(success = FALSE, error = err_msg)
})
',
output_file_path,
email_address,
maximumknotallowed,
permuiter
  )
  
  writeLines(lamian_r_code, r_script_path)
  
  # generate qsub
  log_path <- file.path(output_file_path, "output.log")
  qsub_script_path <- file.path(output_file_path, "lamian.qsub")
  
  qsub_content <- paste0(
    "#!/bin/bash\n",
    "#$ -N lamian\n",
    "#$ -cwd\n",
    "#$ -j y\n",
    "#$ -o ", log_path, "\n",
    "#$ -pe omp 28\n",
    "#$ -l h_rt=24:00:00\n",
    "#$ -l mem_per_core=18G\n",
    "#$ -V\n",
    "#$ -P ", project_name, "\n",
    "\n",
    "module load R/4.4.3\n",
    "source /projectnb/wax-es/Bingtian/Lamian/miniforge_env/etc/profile.d/conda.sh\n",
    "conda activate /projectnb/wax-es/Bingtian/Lamian/miniforge_env\n",
    "\n",
    "Rscript ", r_script_path, "\n"
  )
  writeLines(qsub_content, qsub_script_path)
  
  # submit qsub and return job id
  qsub_output <- system(paste("qsub", qsub_script_path), intern = TRUE)
  job_id <- sub(".*?([0-9]+).*", "\\1", qsub_output)
  return(job_id)
}